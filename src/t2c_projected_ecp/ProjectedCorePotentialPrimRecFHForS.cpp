#include "ProjectedCorePotentialPrimRecFHForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_fh_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_fh_s_0_0_0,
                                        const size_t idx_ph_s_0_0_0,
                                        const size_t idx_dh_s_0_0_0,
                                        const size_t idx_ph_s_1_0_0,
                                        const size_t idx_dh_s_1_0_0,
                                        const int p,
                                        const size_t idx_ph_s_0_0_1,
                                        const size_t idx_dh_s_0_0_1,
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

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0);

    auto tg_x_xxxxy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 1);

    auto tg_x_xxxxz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 2);

    auto tg_x_xxxyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 3);

    auto tg_x_xxxyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 4);

    auto tg_x_xxxzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 5);

    auto tg_x_xxyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 6);

    auto tg_x_xxyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 7);

    auto tg_x_xxyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 8);

    auto tg_x_xxzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 9);

    auto tg_x_xyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 10);

    auto tg_x_xyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 11);

    auto tg_x_xyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 12);

    auto tg_x_xyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 13);

    auto tg_x_xzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 14);

    auto tg_x_yyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 15);

    auto tg_x_yyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 16);

    auto tg_x_yyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 17);

    auto tg_x_yyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 18);

    auto tg_x_yzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 19);

    auto tg_x_zzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 20);

    auto tg_y_xxxxx_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 21);

    auto tg_y_xxxxy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 22);

    auto tg_y_xxxxz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 23);

    auto tg_y_xxxyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 24);

    auto tg_y_xxxyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 25);

    auto tg_y_xxxzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 26);

    auto tg_y_xxyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 27);

    auto tg_y_xxyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 28);

    auto tg_y_xxyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 29);

    auto tg_y_xxzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 30);

    auto tg_y_xyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 31);

    auto tg_y_xyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 32);

    auto tg_y_xyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 33);

    auto tg_y_xyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 34);

    auto tg_y_xzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 35);

    auto tg_y_yyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 36);

    auto tg_y_yyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 37);

    auto tg_y_yyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 38);

    auto tg_y_yyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 39);

    auto tg_y_yzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 40);

    auto tg_y_zzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 41);

    auto tg_z_xxxxx_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 42);

    auto tg_z_xxxxy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 43);

    auto tg_z_xxxxz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 44);

    auto tg_z_xxxyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 45);

    auto tg_z_xxxyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 46);

    auto tg_z_xxxzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 47);

    auto tg_z_xxyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 48);

    auto tg_z_xxyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 49);

    auto tg_z_xxyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 50);

    auto tg_z_xxzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 51);

    auto tg_z_xyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 52);

    auto tg_z_xyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 53);

    auto tg_z_xyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 54);

    auto tg_z_xyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 55);

    auto tg_z_xzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 56);

    auto tg_z_yyyyy_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 57);

    auto tg_z_yyyyz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 58);

    auto tg_z_yyyzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 59);

    auto tg_z_yyzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 60);

    auto tg_z_yzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 61);

    auto tg_z_zzzzz_s_0_0_0 = pbuffer.data(idx_ph_s_0_0_0 + 62);

    // Set up components of auxiliary buffer : DH

    auto tg_xx_xxxxx_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0);

    auto tg_xx_xxxxy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 1);

    auto tg_xx_xxxxz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 2);

    auto tg_xx_xxxyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 3);

    auto tg_xx_xxxyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 4);

    auto tg_xx_xxxzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 5);

    auto tg_xx_xxyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 6);

    auto tg_xx_xxyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 7);

    auto tg_xx_xxyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 8);

    auto tg_xx_xxzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 9);

    auto tg_xx_xyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 10);

    auto tg_xx_xyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 11);

    auto tg_xx_xyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 12);

    auto tg_xx_xyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 13);

    auto tg_xx_xzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 14);

    auto tg_xx_yyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 15);

    auto tg_xx_yyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 16);

    auto tg_xx_yyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 17);

    auto tg_xx_yyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 18);

    auto tg_xx_yzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 19);

    auto tg_xx_zzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 20);

    auto tg_xy_xxxxx_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 21);

    auto tg_xy_xxxxy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 22);

    auto tg_xy_xxxxz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 23);

    auto tg_xy_xxxyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 24);

    auto tg_xy_xxxyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 25);

    auto tg_xy_xxxzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 26);

    auto tg_xy_xxyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 27);

    auto tg_xy_xxyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 28);

    auto tg_xy_xxyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 29);

    auto tg_xy_xxzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 30);

    auto tg_xy_xyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 31);

    auto tg_xy_xyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 32);

    auto tg_xy_xyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 33);

    auto tg_xy_xyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 34);

    auto tg_xy_xzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 35);

    auto tg_xy_yyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 36);

    auto tg_xy_yyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 37);

    auto tg_xy_yyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 38);

    auto tg_xy_yyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 39);

    auto tg_xy_yzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 40);

    auto tg_xy_zzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 41);

    auto tg_xz_xxxxx_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 42);

    auto tg_xz_xxxxy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 43);

    auto tg_xz_xxxxz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 44);

    auto tg_xz_xxxyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 45);

    auto tg_xz_xxxyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 46);

    auto tg_xz_xxxzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 47);

    auto tg_xz_xxyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 48);

    auto tg_xz_xxyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 49);

    auto tg_xz_xxyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 50);

    auto tg_xz_xxzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 51);

    auto tg_xz_xyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 52);

    auto tg_xz_xyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 53);

    auto tg_xz_xyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 54);

    auto tg_xz_xyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 55);

    auto tg_xz_xzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 56);

    auto tg_xz_yyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 57);

    auto tg_xz_yyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 58);

    auto tg_xz_yyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 59);

    auto tg_xz_yyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 60);

    auto tg_xz_yzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 61);

    auto tg_xz_zzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 62);

    auto tg_yy_xxxxx_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 63);

    auto tg_yy_xxxxy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 64);

    auto tg_yy_xxxxz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 65);

    auto tg_yy_xxxyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 66);

    auto tg_yy_xxxyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 67);

    auto tg_yy_xxxzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 68);

    auto tg_yy_xxyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 69);

    auto tg_yy_xxyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 70);

    auto tg_yy_xxyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 71);

    auto tg_yy_xxzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 72);

    auto tg_yy_xyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 73);

    auto tg_yy_xyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 74);

    auto tg_yy_xyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 75);

    auto tg_yy_xyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 76);

    auto tg_yy_xzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 77);

    auto tg_yy_yyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 78);

    auto tg_yy_yyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 79);

    auto tg_yy_yyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 80);

    auto tg_yy_yyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 81);

    auto tg_yy_yzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 82);

    auto tg_yy_zzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 83);

    auto tg_yz_xxxxx_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 84);

    auto tg_yz_xxxxy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 85);

    auto tg_yz_xxxxz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 86);

    auto tg_yz_xxxyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 87);

    auto tg_yz_xxxyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 88);

    auto tg_yz_xxxzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 89);

    auto tg_yz_xxyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 90);

    auto tg_yz_xxyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 91);

    auto tg_yz_xxyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 92);

    auto tg_yz_xxzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 93);

    auto tg_yz_xyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 94);

    auto tg_yz_xyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 95);

    auto tg_yz_xyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 96);

    auto tg_yz_xyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 97);

    auto tg_yz_xzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 98);

    auto tg_yz_yyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 99);

    auto tg_yz_yyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 100);

    auto tg_yz_yyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 101);

    auto tg_yz_yyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 102);

    auto tg_yz_yzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 103);

    auto tg_yz_zzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 104);

    auto tg_zz_xxxxx_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 105);

    auto tg_zz_xxxxy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 106);

    auto tg_zz_xxxxz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 107);

    auto tg_zz_xxxyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 108);

    auto tg_zz_xxxyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 109);

    auto tg_zz_xxxzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 110);

    auto tg_zz_xxyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 111);

    auto tg_zz_xxyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 112);

    auto tg_zz_xxyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 113);

    auto tg_zz_xxzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 114);

    auto tg_zz_xyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 115);

    auto tg_zz_xyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 116);

    auto tg_zz_xyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 117);

    auto tg_zz_xyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 118);

    auto tg_zz_xzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 119);

    auto tg_zz_yyyyy_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 120);

    auto tg_zz_yyyyz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 121);

    auto tg_zz_yyyzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 122);

    auto tg_zz_yyzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 123);

    auto tg_zz_yzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 124);

    auto tg_zz_zzzzz_s_0_0_0 = pbuffer.data(idx_dh_s_0_0_0 + 125);

    // Set up components of auxiliary buffer : PH

    auto tg_x_xxxxx_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0);

    auto tg_x_xxxxy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 1);

    auto tg_x_xxxxz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 2);

    auto tg_x_xxxyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 3);

    auto tg_x_xxxyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 4);

    auto tg_x_xxxzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 5);

    auto tg_x_xxyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 6);

    auto tg_x_xxyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 7);

    auto tg_x_xxyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 8);

    auto tg_x_xxzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 9);

    auto tg_x_xyyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 10);

    auto tg_x_xyyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 11);

    auto tg_x_xyyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 12);

    auto tg_x_xyzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 13);

    auto tg_x_xzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 14);

    auto tg_x_yyyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 15);

    auto tg_x_yyyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 16);

    auto tg_x_yyyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 17);

    auto tg_x_yyzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 18);

    auto tg_x_yzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 19);

    auto tg_x_zzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 20);

    auto tg_y_xxxxx_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 21);

    auto tg_y_xxxxy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 22);

    auto tg_y_xxxxz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 23);

    auto tg_y_xxxyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 24);

    auto tg_y_xxxyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 25);

    auto tg_y_xxxzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 26);

    auto tg_y_xxyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 27);

    auto tg_y_xxyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 28);

    auto tg_y_xxyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 29);

    auto tg_y_xxzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 30);

    auto tg_y_xyyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 31);

    auto tg_y_xyyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 32);

    auto tg_y_xyyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 33);

    auto tg_y_xyzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 34);

    auto tg_y_xzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 35);

    auto tg_y_yyyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 36);

    auto tg_y_yyyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 37);

    auto tg_y_yyyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 38);

    auto tg_y_yyzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 39);

    auto tg_y_yzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 40);

    auto tg_y_zzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 41);

    auto tg_z_xxxxx_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 42);

    auto tg_z_xxxxy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 43);

    auto tg_z_xxxxz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 44);

    auto tg_z_xxxyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 45);

    auto tg_z_xxxyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 46);

    auto tg_z_xxxzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 47);

    auto tg_z_xxyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 48);

    auto tg_z_xxyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 49);

    auto tg_z_xxyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 50);

    auto tg_z_xxzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 51);

    auto tg_z_xyyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 52);

    auto tg_z_xyyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 53);

    auto tg_z_xyyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 54);

    auto tg_z_xyzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 55);

    auto tg_z_xzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 56);

    auto tg_z_yyyyy_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 57);

    auto tg_z_yyyyz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 58);

    auto tg_z_yyyzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 59);

    auto tg_z_yyzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 60);

    auto tg_z_yzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 61);

    auto tg_z_zzzzz_s_1_0_0 = pbuffer.data(idx_ph_s_1_0_0 + 62);

    // Set up components of auxiliary buffer : DH

    auto tg_xx_xxxxx_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0);

    auto tg_xx_xxxxy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 1);

    auto tg_xx_xxxxz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 2);

    auto tg_xx_xxxyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 3);

    auto tg_xx_xxxyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 4);

    auto tg_xx_xxxzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 5);

    auto tg_xx_xxyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 6);

    auto tg_xx_xxyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 7);

    auto tg_xx_xxyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 8);

    auto tg_xx_xxzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 9);

    auto tg_xx_xyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 10);

    auto tg_xx_xyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 11);

    auto tg_xx_xyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 12);

    auto tg_xx_xyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 13);

    auto tg_xx_xzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 14);

    auto tg_xx_yyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 15);

    auto tg_xx_yyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 16);

    auto tg_xx_yyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 17);

    auto tg_xx_yyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 18);

    auto tg_xx_yzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 19);

    auto tg_xx_zzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 20);

    auto tg_xy_xxxxx_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 21);

    auto tg_xy_xxxxy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 22);

    auto tg_xy_xxxxz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 23);

    auto tg_xy_xxxyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 24);

    auto tg_xy_xxxyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 25);

    auto tg_xy_xxxzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 26);

    auto tg_xy_xxyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 27);

    auto tg_xy_xxyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 28);

    auto tg_xy_xxyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 29);

    auto tg_xy_xxzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 30);

    auto tg_xy_xyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 31);

    auto tg_xy_xyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 32);

    auto tg_xy_xyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 33);

    auto tg_xy_xyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 34);

    auto tg_xy_xzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 35);

    auto tg_xy_yyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 36);

    auto tg_xy_yyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 37);

    auto tg_xy_yyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 38);

    auto tg_xy_yyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 39);

    auto tg_xy_yzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 40);

    auto tg_xy_zzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 41);

    auto tg_xz_xxxxx_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 42);

    auto tg_xz_xxxxy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 43);

    auto tg_xz_xxxxz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 44);

    auto tg_xz_xxxyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 45);

    auto tg_xz_xxxyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 46);

    auto tg_xz_xxxzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 47);

    auto tg_xz_xxyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 48);

    auto tg_xz_xxyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 49);

    auto tg_xz_xxyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 50);

    auto tg_xz_xxzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 51);

    auto tg_xz_xyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 52);

    auto tg_xz_xyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 53);

    auto tg_xz_xyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 54);

    auto tg_xz_xyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 55);

    auto tg_xz_xzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 56);

    auto tg_xz_yyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 57);

    auto tg_xz_yyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 58);

    auto tg_xz_yyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 59);

    auto tg_xz_yyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 60);

    auto tg_xz_yzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 61);

    auto tg_xz_zzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 62);

    auto tg_yy_xxxxx_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 63);

    auto tg_yy_xxxxy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 64);

    auto tg_yy_xxxxz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 65);

    auto tg_yy_xxxyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 66);

    auto tg_yy_xxxyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 67);

    auto tg_yy_xxxzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 68);

    auto tg_yy_xxyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 69);

    auto tg_yy_xxyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 70);

    auto tg_yy_xxyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 71);

    auto tg_yy_xxzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 72);

    auto tg_yy_xyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 73);

    auto tg_yy_xyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 74);

    auto tg_yy_xyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 75);

    auto tg_yy_xyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 76);

    auto tg_yy_xzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 77);

    auto tg_yy_yyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 78);

    auto tg_yy_yyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 79);

    auto tg_yy_yyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 80);

    auto tg_yy_yyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 81);

    auto tg_yy_yzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 82);

    auto tg_yy_zzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 83);

    auto tg_yz_xxxxx_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 84);

    auto tg_yz_xxxxy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 85);

    auto tg_yz_xxxxz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 86);

    auto tg_yz_xxxyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 87);

    auto tg_yz_xxxyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 88);

    auto tg_yz_xxxzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 89);

    auto tg_yz_xxyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 90);

    auto tg_yz_xxyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 91);

    auto tg_yz_xxyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 92);

    auto tg_yz_xxzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 93);

    auto tg_yz_xyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 94);

    auto tg_yz_xyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 95);

    auto tg_yz_xyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 96);

    auto tg_yz_xyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 97);

    auto tg_yz_xzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 98);

    auto tg_yz_yyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 99);

    auto tg_yz_yyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 100);

    auto tg_yz_yyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 101);

    auto tg_yz_yyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 102);

    auto tg_yz_yzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 103);

    auto tg_yz_zzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 104);

    auto tg_zz_xxxxx_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 105);

    auto tg_zz_xxxxy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 106);

    auto tg_zz_xxxxz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 107);

    auto tg_zz_xxxyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 108);

    auto tg_zz_xxxyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 109);

    auto tg_zz_xxxzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 110);

    auto tg_zz_xxyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 111);

    auto tg_zz_xxyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 112);

    auto tg_zz_xxyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 113);

    auto tg_zz_xxzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 114);

    auto tg_zz_xyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 115);

    auto tg_zz_xyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 116);

    auto tg_zz_xyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 117);

    auto tg_zz_xyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 118);

    auto tg_zz_xzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 119);

    auto tg_zz_yyyyy_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 120);

    auto tg_zz_yyyyz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 121);

    auto tg_zz_yyyzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 122);

    auto tg_zz_yyzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 123);

    auto tg_zz_yzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 124);

    auto tg_zz_zzzzz_s_1_0_0 = pbuffer.data(idx_dh_s_1_0_0 + 125);

    // Set up components of targeted buffer : FH

    auto tg_xxx_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0);

    auto tg_xxx_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 1);

    auto tg_xxx_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 2);

    auto tg_xxx_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 3);

    auto tg_xxx_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 4);

    auto tg_xxx_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 5);

    auto tg_xxx_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 6);

    auto tg_xxx_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 7);

    auto tg_xxx_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 8);

    auto tg_xxx_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 9);

    auto tg_xxx_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 10);

    auto tg_xxx_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 11);

    auto tg_xxx_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 12);

    auto tg_xxx_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 13);

    auto tg_xxx_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 14);

    auto tg_xxx_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 15);

    auto tg_xxx_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 16);

    auto tg_xxx_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 17);

    auto tg_xxx_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 18);

    auto tg_xxx_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 19);

    auto tg_xxx_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 20);

    auto tg_xxy_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 21);

    auto tg_xxy_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 22);

    auto tg_xxy_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 23);

    auto tg_xxy_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 24);

    auto tg_xxy_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 25);

    auto tg_xxy_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 26);

    auto tg_xxy_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 27);

    auto tg_xxy_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 28);

    auto tg_xxy_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 29);

    auto tg_xxy_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 30);

    auto tg_xxy_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 31);

    auto tg_xxy_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 32);

    auto tg_xxy_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 33);

    auto tg_xxy_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 34);

    auto tg_xxy_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 35);

    auto tg_xxy_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 36);

    auto tg_xxy_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 37);

    auto tg_xxy_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 38);

    auto tg_xxy_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 39);

    auto tg_xxy_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 40);

    auto tg_xxy_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 41);

    auto tg_xxz_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 42);

    auto tg_xxz_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 43);

    auto tg_xxz_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 44);

    auto tg_xxz_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 45);

    auto tg_xxz_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 46);

    auto tg_xxz_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 47);

    auto tg_xxz_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 48);

    auto tg_xxz_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 49);

    auto tg_xxz_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 50);

    auto tg_xxz_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 51);

    auto tg_xxz_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 52);

    auto tg_xxz_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 53);

    auto tg_xxz_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 54);

    auto tg_xxz_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 55);

    auto tg_xxz_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 56);

    auto tg_xxz_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 57);

    auto tg_xxz_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 58);

    auto tg_xxz_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 59);

    auto tg_xxz_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 60);

    auto tg_xxz_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 61);

    auto tg_xxz_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 62);

    auto tg_xyy_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 63);

    auto tg_xyy_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 64);

    auto tg_xyy_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 65);

    auto tg_xyy_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 66);

    auto tg_xyy_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 67);

    auto tg_xyy_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 68);

    auto tg_xyy_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 69);

    auto tg_xyy_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 70);

    auto tg_xyy_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 71);

    auto tg_xyy_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 72);

    auto tg_xyy_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 73);

    auto tg_xyy_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 74);

    auto tg_xyy_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 75);

    auto tg_xyy_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 76);

    auto tg_xyy_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 77);

    auto tg_xyy_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 78);

    auto tg_xyy_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 79);

    auto tg_xyy_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 80);

    auto tg_xyy_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 81);

    auto tg_xyy_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 82);

    auto tg_xyy_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 83);

    auto tg_xyz_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 84);

    auto tg_xyz_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 85);

    auto tg_xyz_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 86);

    auto tg_xyz_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 87);

    auto tg_xyz_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 88);

    auto tg_xyz_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 89);

    auto tg_xyz_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 90);

    auto tg_xyz_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 91);

    auto tg_xyz_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 92);

    auto tg_xyz_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 93);

    auto tg_xyz_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 94);

    auto tg_xyz_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 95);

    auto tg_xyz_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 96);

    auto tg_xyz_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 97);

    auto tg_xyz_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 98);

    auto tg_xyz_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 99);

    auto tg_xyz_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 100);

    auto tg_xyz_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 101);

    auto tg_xyz_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 102);

    auto tg_xyz_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 103);

    auto tg_xyz_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 104);

    auto tg_xzz_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 105);

    auto tg_xzz_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 106);

    auto tg_xzz_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 107);

    auto tg_xzz_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 108);

    auto tg_xzz_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 109);

    auto tg_xzz_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 110);

    auto tg_xzz_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 111);

    auto tg_xzz_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 112);

    auto tg_xzz_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 113);

    auto tg_xzz_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 114);

    auto tg_xzz_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 115);

    auto tg_xzz_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 116);

    auto tg_xzz_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 117);

    auto tg_xzz_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 118);

    auto tg_xzz_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 119);

    auto tg_xzz_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 120);

    auto tg_xzz_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 121);

    auto tg_xzz_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 122);

    auto tg_xzz_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 123);

    auto tg_xzz_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 124);

    auto tg_xzz_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 125);

    auto tg_yyy_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 126);

    auto tg_yyy_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 127);

    auto tg_yyy_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 128);

    auto tg_yyy_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 129);

    auto tg_yyy_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 130);

    auto tg_yyy_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 131);

    auto tg_yyy_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 132);

    auto tg_yyy_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 133);

    auto tg_yyy_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 134);

    auto tg_yyy_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 135);

    auto tg_yyy_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 136);

    auto tg_yyy_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 137);

    auto tg_yyy_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 138);

    auto tg_yyy_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 139);

    auto tg_yyy_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 140);

    auto tg_yyy_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 141);

    auto tg_yyy_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 142);

    auto tg_yyy_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 143);

    auto tg_yyy_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 144);

    auto tg_yyy_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 145);

    auto tg_yyy_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 146);

    auto tg_yyz_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 147);

    auto tg_yyz_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 148);

    auto tg_yyz_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 149);

    auto tg_yyz_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 150);

    auto tg_yyz_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 151);

    auto tg_yyz_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 152);

    auto tg_yyz_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 153);

    auto tg_yyz_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 154);

    auto tg_yyz_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 155);

    auto tg_yyz_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 156);

    auto tg_yyz_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 157);

    auto tg_yyz_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 158);

    auto tg_yyz_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 159);

    auto tg_yyz_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 160);

    auto tg_yyz_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 161);

    auto tg_yyz_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 162);

    auto tg_yyz_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 163);

    auto tg_yyz_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 164);

    auto tg_yyz_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 165);

    auto tg_yyz_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 166);

    auto tg_yyz_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 167);

    auto tg_yzz_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 168);

    auto tg_yzz_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 169);

    auto tg_yzz_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 170);

    auto tg_yzz_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 171);

    auto tg_yzz_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 172);

    auto tg_yzz_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 173);

    auto tg_yzz_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 174);

    auto tg_yzz_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 175);

    auto tg_yzz_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 176);

    auto tg_yzz_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 177);

    auto tg_yzz_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 178);

    auto tg_yzz_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 179);

    auto tg_yzz_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 180);

    auto tg_yzz_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 181);

    auto tg_yzz_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 182);

    auto tg_yzz_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 183);

    auto tg_yzz_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 184);

    auto tg_yzz_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 185);

    auto tg_yzz_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 186);

    auto tg_yzz_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 187);

    auto tg_yzz_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 188);

    auto tg_zzz_xxxxx_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 189);

    auto tg_zzz_xxxxy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 190);

    auto tg_zzz_xxxxz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 191);

    auto tg_zzz_xxxyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 192);

    auto tg_zzz_xxxyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 193);

    auto tg_zzz_xxxzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 194);

    auto tg_zzz_xxyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 195);

    auto tg_zzz_xxyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 196);

    auto tg_zzz_xxyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 197);

    auto tg_zzz_xxzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 198);

    auto tg_zzz_xyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 199);

    auto tg_zzz_xyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 200);

    auto tg_zzz_xyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 201);

    auto tg_zzz_xyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 202);

    auto tg_zzz_xzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 203);

    auto tg_zzz_yyyyy_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 204);

    auto tg_zzz_yyyyz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 205);

    auto tg_zzz_yyyzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 206);

    auto tg_zzz_yyzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 207);

    auto tg_zzz_yzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 208);

    auto tg_zzz_zzzzz_s_0_0_0 = pbuffer.data(idx_fh_s_0_0_0 + 209);

    #pragma omp simd aligned(b_exps, tg_x_xxxxx_s_0_0_0, tg_x_xxxxx_s_1_0_0, tg_x_xxxxy_s_0_0_0, tg_x_xxxxy_s_1_0_0, tg_x_xxxxz_s_0_0_0, tg_x_xxxxz_s_1_0_0, tg_x_xxxyy_s_0_0_0, tg_x_xxxyy_s_1_0_0, tg_x_xxxyz_s_0_0_0, tg_x_xxxyz_s_1_0_0, tg_x_xxxzz_s_0_0_0, tg_x_xxxzz_s_1_0_0, tg_x_xxyyy_s_0_0_0, tg_x_xxyyy_s_1_0_0, tg_x_xxyyz_s_0_0_0, tg_x_xxyyz_s_1_0_0, tg_x_xxyzz_s_0_0_0, tg_x_xxyzz_s_1_0_0, tg_x_xxzzz_s_0_0_0, tg_x_xxzzz_s_1_0_0, tg_x_xyyyy_s_0_0_0, tg_x_xyyyy_s_1_0_0, tg_x_xyyyz_s_0_0_0, tg_x_xyyyz_s_1_0_0, tg_x_xyyzz_s_0_0_0, tg_x_xyyzz_s_1_0_0, tg_x_xyzzz_s_0_0_0, tg_x_xyzzz_s_1_0_0, tg_x_xzzzz_s_0_0_0, tg_x_xzzzz_s_1_0_0, tg_x_yyyyy_s_0_0_0, tg_x_yyyyy_s_1_0_0, tg_x_yyyyz_s_0_0_0, tg_x_yyyyz_s_1_0_0, tg_x_yyyzz_s_0_0_0, tg_x_yyyzz_s_1_0_0, tg_x_yyzzz_s_0_0_0, tg_x_yyzzz_s_1_0_0, tg_x_yzzzz_s_0_0_0, tg_x_yzzzz_s_1_0_0, tg_x_zzzzz_s_0_0_0, tg_x_zzzzz_s_1_0_0, tg_xx_xxxxx_s_0_0_0, tg_xx_xxxxx_s_1_0_0, tg_xx_xxxxy_s_0_0_0, tg_xx_xxxxy_s_1_0_0, tg_xx_xxxxz_s_0_0_0, tg_xx_xxxxz_s_1_0_0, tg_xx_xxxyy_s_0_0_0, tg_xx_xxxyy_s_1_0_0, tg_xx_xxxyz_s_0_0_0, tg_xx_xxxyz_s_1_0_0, tg_xx_xxxzz_s_0_0_0, tg_xx_xxxzz_s_1_0_0, tg_xx_xxyyy_s_0_0_0, tg_xx_xxyyy_s_1_0_0, tg_xx_xxyyz_s_0_0_0, tg_xx_xxyyz_s_1_0_0, tg_xx_xxyzz_s_0_0_0, tg_xx_xxyzz_s_1_0_0, tg_xx_xxzzz_s_0_0_0, tg_xx_xxzzz_s_1_0_0, tg_xx_xyyyy_s_0_0_0, tg_xx_xyyyy_s_1_0_0, tg_xx_xyyyz_s_0_0_0, tg_xx_xyyyz_s_1_0_0, tg_xx_xyyzz_s_0_0_0, tg_xx_xyyzz_s_1_0_0, tg_xx_xyzzz_s_0_0_0, tg_xx_xyzzz_s_1_0_0, tg_xx_xzzzz_s_0_0_0, tg_xx_xzzzz_s_1_0_0, tg_xx_yyyyy_s_0_0_0, tg_xx_yyyyy_s_1_0_0, tg_xx_yyyyz_s_0_0_0, tg_xx_yyyyz_s_1_0_0, tg_xx_yyyzz_s_0_0_0, tg_xx_yyyzz_s_1_0_0, tg_xx_yyzzz_s_0_0_0, tg_xx_yyzzz_s_1_0_0, tg_xx_yzzzz_s_0_0_0, tg_xx_yzzzz_s_1_0_0, tg_xx_zzzzz_s_0_0_0, tg_xx_zzzzz_s_1_0_0, tg_xxx_xxxxx_s_0_0_0, tg_xxx_xxxxy_s_0_0_0, tg_xxx_xxxxz_s_0_0_0, tg_xxx_xxxyy_s_0_0_0, tg_xxx_xxxyz_s_0_0_0, tg_xxx_xxxzz_s_0_0_0, tg_xxx_xxyyy_s_0_0_0, tg_xxx_xxyyz_s_0_0_0, tg_xxx_xxyzz_s_0_0_0, tg_xxx_xxzzz_s_0_0_0, tg_xxx_xyyyy_s_0_0_0, tg_xxx_xyyyz_s_0_0_0, tg_xxx_xyyzz_s_0_0_0, tg_xxx_xyzzz_s_0_0_0, tg_xxx_xzzzz_s_0_0_0, tg_xxx_yyyyy_s_0_0_0, tg_xxx_yyyyz_s_0_0_0, tg_xxx_yyyzz_s_0_0_0, tg_xxx_yyzzz_s_0_0_0, tg_xxx_yzzzz_s_0_0_0, tg_xxx_zzzzz_s_0_0_0, tg_xxy_xxxxx_s_0_0_0, tg_xxy_xxxxy_s_0_0_0, tg_xxy_xxxxz_s_0_0_0, tg_xxy_xxxyy_s_0_0_0, tg_xxy_xxxyz_s_0_0_0, tg_xxy_xxxzz_s_0_0_0, tg_xxy_xxyyy_s_0_0_0, tg_xxy_xxyyz_s_0_0_0, tg_xxy_xxyzz_s_0_0_0, tg_xxy_xxzzz_s_0_0_0, tg_xxy_xyyyy_s_0_0_0, tg_xxy_xyyyz_s_0_0_0, tg_xxy_xyyzz_s_0_0_0, tg_xxy_xyzzz_s_0_0_0, tg_xxy_xzzzz_s_0_0_0, tg_xxy_yyyyy_s_0_0_0, tg_xxy_yyyyz_s_0_0_0, tg_xxy_yyyzz_s_0_0_0, tg_xxy_yyzzz_s_0_0_0, tg_xxy_yzzzz_s_0_0_0, tg_xxy_zzzzz_s_0_0_0, tg_xxz_xxxxx_s_0_0_0, tg_xxz_xxxxy_s_0_0_0, tg_xxz_xxxxz_s_0_0_0, tg_xxz_xxxyy_s_0_0_0, tg_xxz_xxxyz_s_0_0_0, tg_xxz_xxxzz_s_0_0_0, tg_xxz_xxyyy_s_0_0_0, tg_xxz_xxyyz_s_0_0_0, tg_xxz_xxyzz_s_0_0_0, tg_xxz_xxzzz_s_0_0_0, tg_xxz_xyyyy_s_0_0_0, tg_xxz_xyyyz_s_0_0_0, tg_xxz_xyyzz_s_0_0_0, tg_xxz_xyzzz_s_0_0_0, tg_xxz_xzzzz_s_0_0_0, tg_xxz_yyyyy_s_0_0_0, tg_xxz_yyyyz_s_0_0_0, tg_xxz_yyyzz_s_0_0_0, tg_xxz_yyzzz_s_0_0_0, tg_xxz_yzzzz_s_0_0_0, tg_xxz_zzzzz_s_0_0_0, tg_xyy_xxxxx_s_0_0_0, tg_xyy_xxxxy_s_0_0_0, tg_xyy_xxxxz_s_0_0_0, tg_xyy_xxxyy_s_0_0_0, tg_xyy_xxxyz_s_0_0_0, tg_xyy_xxxzz_s_0_0_0, tg_xyy_xxyyy_s_0_0_0, tg_xyy_xxyyz_s_0_0_0, tg_xyy_xxyzz_s_0_0_0, tg_xyy_xxzzz_s_0_0_0, tg_xyy_xyyyy_s_0_0_0, tg_xyy_xyyyz_s_0_0_0, tg_xyy_xyyzz_s_0_0_0, tg_xyy_xyzzz_s_0_0_0, tg_xyy_xzzzz_s_0_0_0, tg_xyy_yyyyy_s_0_0_0, tg_xyy_yyyyz_s_0_0_0, tg_xyy_yyyzz_s_0_0_0, tg_xyy_yyzzz_s_0_0_0, tg_xyy_yzzzz_s_0_0_0, tg_xyy_zzzzz_s_0_0_0, tg_xyz_xxxxx_s_0_0_0, tg_xyz_xxxxy_s_0_0_0, tg_xyz_xxxxz_s_0_0_0, tg_xyz_xxxyy_s_0_0_0, tg_xyz_xxxyz_s_0_0_0, tg_xyz_xxxzz_s_0_0_0, tg_xyz_xxyyy_s_0_0_0, tg_xyz_xxyyz_s_0_0_0, tg_xyz_xxyzz_s_0_0_0, tg_xyz_xxzzz_s_0_0_0, tg_xyz_xyyyy_s_0_0_0, tg_xyz_xyyyz_s_0_0_0, tg_xyz_xyyzz_s_0_0_0, tg_xyz_xyzzz_s_0_0_0, tg_xyz_xzzzz_s_0_0_0, tg_xyz_yyyyy_s_0_0_0, tg_xyz_yyyyz_s_0_0_0, tg_xyz_yyyzz_s_0_0_0, tg_xyz_yyzzz_s_0_0_0, tg_xyz_yzzzz_s_0_0_0, tg_xyz_zzzzz_s_0_0_0, tg_xzz_xxxxx_s_0_0_0, tg_xzz_xxxxy_s_0_0_0, tg_xzz_xxxxz_s_0_0_0, tg_xzz_xxxyy_s_0_0_0, tg_xzz_xxxyz_s_0_0_0, tg_xzz_xxxzz_s_0_0_0, tg_xzz_xxyyy_s_0_0_0, tg_xzz_xxyyz_s_0_0_0, tg_xzz_xxyzz_s_0_0_0, tg_xzz_xxzzz_s_0_0_0, tg_xzz_xyyyy_s_0_0_0, tg_xzz_xyyyz_s_0_0_0, tg_xzz_xyyzz_s_0_0_0, tg_xzz_xyzzz_s_0_0_0, tg_xzz_xzzzz_s_0_0_0, tg_xzz_yyyyy_s_0_0_0, tg_xzz_yyyyz_s_0_0_0, tg_xzz_yyyzz_s_0_0_0, tg_xzz_yyzzz_s_0_0_0, tg_xzz_yzzzz_s_0_0_0, tg_xzz_zzzzz_s_0_0_0, tg_y_xxxxx_s_0_0_0, tg_y_xxxxx_s_1_0_0, tg_y_xxxxy_s_0_0_0, tg_y_xxxxy_s_1_0_0, tg_y_xxxxz_s_0_0_0, tg_y_xxxxz_s_1_0_0, tg_y_xxxyy_s_0_0_0, tg_y_xxxyy_s_1_0_0, tg_y_xxxyz_s_0_0_0, tg_y_xxxyz_s_1_0_0, tg_y_xxxzz_s_0_0_0, tg_y_xxxzz_s_1_0_0, tg_y_xxyyy_s_0_0_0, tg_y_xxyyy_s_1_0_0, tg_y_xxyyz_s_0_0_0, tg_y_xxyyz_s_1_0_0, tg_y_xxyzz_s_0_0_0, tg_y_xxyzz_s_1_0_0, tg_y_xxzzz_s_0_0_0, tg_y_xxzzz_s_1_0_0, tg_y_xyyyy_s_0_0_0, tg_y_xyyyy_s_1_0_0, tg_y_xyyyz_s_0_0_0, tg_y_xyyyz_s_1_0_0, tg_y_xyyzz_s_0_0_0, tg_y_xyyzz_s_1_0_0, tg_y_xyzzz_s_0_0_0, tg_y_xyzzz_s_1_0_0, tg_y_xzzzz_s_0_0_0, tg_y_xzzzz_s_1_0_0, tg_y_yyyyy_s_0_0_0, tg_y_yyyyy_s_1_0_0, tg_y_yyyyz_s_0_0_0, tg_y_yyyyz_s_1_0_0, tg_y_yyyzz_s_0_0_0, tg_y_yyyzz_s_1_0_0, tg_y_yyzzz_s_0_0_0, tg_y_yyzzz_s_1_0_0, tg_y_yzzzz_s_0_0_0, tg_y_yzzzz_s_1_0_0, tg_y_zzzzz_s_0_0_0, tg_y_zzzzz_s_1_0_0, tg_yy_xxxxx_s_0_0_0, tg_yy_xxxxx_s_1_0_0, tg_yy_xxxxy_s_0_0_0, tg_yy_xxxxy_s_1_0_0, tg_yy_xxxxz_s_0_0_0, tg_yy_xxxxz_s_1_0_0, tg_yy_xxxyy_s_0_0_0, tg_yy_xxxyy_s_1_0_0, tg_yy_xxxyz_s_0_0_0, tg_yy_xxxyz_s_1_0_0, tg_yy_xxxzz_s_0_0_0, tg_yy_xxxzz_s_1_0_0, tg_yy_xxyyy_s_0_0_0, tg_yy_xxyyy_s_1_0_0, tg_yy_xxyyz_s_0_0_0, tg_yy_xxyyz_s_1_0_0, tg_yy_xxyzz_s_0_0_0, tg_yy_xxyzz_s_1_0_0, tg_yy_xxzzz_s_0_0_0, tg_yy_xxzzz_s_1_0_0, tg_yy_xyyyy_s_0_0_0, tg_yy_xyyyy_s_1_0_0, tg_yy_xyyyz_s_0_0_0, tg_yy_xyyyz_s_1_0_0, tg_yy_xyyzz_s_0_0_0, tg_yy_xyyzz_s_1_0_0, tg_yy_xyzzz_s_0_0_0, tg_yy_xyzzz_s_1_0_0, tg_yy_xzzzz_s_0_0_0, tg_yy_xzzzz_s_1_0_0, tg_yy_yyyyy_s_0_0_0, tg_yy_yyyyy_s_1_0_0, tg_yy_yyyyz_s_0_0_0, tg_yy_yyyyz_s_1_0_0, tg_yy_yyyzz_s_0_0_0, tg_yy_yyyzz_s_1_0_0, tg_yy_yyzzz_s_0_0_0, tg_yy_yyzzz_s_1_0_0, tg_yy_yzzzz_s_0_0_0, tg_yy_yzzzz_s_1_0_0, tg_yy_zzzzz_s_0_0_0, tg_yy_zzzzz_s_1_0_0, tg_yyy_xxxxx_s_0_0_0, tg_yyy_xxxxy_s_0_0_0, tg_yyy_xxxxz_s_0_0_0, tg_yyy_xxxyy_s_0_0_0, tg_yyy_xxxyz_s_0_0_0, tg_yyy_xxxzz_s_0_0_0, tg_yyy_xxyyy_s_0_0_0, tg_yyy_xxyyz_s_0_0_0, tg_yyy_xxyzz_s_0_0_0, tg_yyy_xxzzz_s_0_0_0, tg_yyy_xyyyy_s_0_0_0, tg_yyy_xyyyz_s_0_0_0, tg_yyy_xyyzz_s_0_0_0, tg_yyy_xyzzz_s_0_0_0, tg_yyy_xzzzz_s_0_0_0, tg_yyy_yyyyy_s_0_0_0, tg_yyy_yyyyz_s_0_0_0, tg_yyy_yyyzz_s_0_0_0, tg_yyy_yyzzz_s_0_0_0, tg_yyy_yzzzz_s_0_0_0, tg_yyy_zzzzz_s_0_0_0, tg_yyz_xxxxx_s_0_0_0, tg_yyz_xxxxy_s_0_0_0, tg_yyz_xxxxz_s_0_0_0, tg_yyz_xxxyy_s_0_0_0, tg_yyz_xxxyz_s_0_0_0, tg_yyz_xxxzz_s_0_0_0, tg_yyz_xxyyy_s_0_0_0, tg_yyz_xxyyz_s_0_0_0, tg_yyz_xxyzz_s_0_0_0, tg_yyz_xxzzz_s_0_0_0, tg_yyz_xyyyy_s_0_0_0, tg_yyz_xyyyz_s_0_0_0, tg_yyz_xyyzz_s_0_0_0, tg_yyz_xyzzz_s_0_0_0, tg_yyz_xzzzz_s_0_0_0, tg_yyz_yyyyy_s_0_0_0, tg_yyz_yyyyz_s_0_0_0, tg_yyz_yyyzz_s_0_0_0, tg_yyz_yyzzz_s_0_0_0, tg_yyz_yzzzz_s_0_0_0, tg_yyz_zzzzz_s_0_0_0, tg_yz_xxxxx_s_0_0_0, tg_yz_xxxxx_s_1_0_0, tg_yz_xxxxy_s_0_0_0, tg_yz_xxxxy_s_1_0_0, tg_yz_xxxxz_s_0_0_0, tg_yz_xxxxz_s_1_0_0, tg_yz_xxxyy_s_0_0_0, tg_yz_xxxyy_s_1_0_0, tg_yz_xxxyz_s_0_0_0, tg_yz_xxxyz_s_1_0_0, tg_yz_xxxzz_s_0_0_0, tg_yz_xxxzz_s_1_0_0, tg_yz_xxyyy_s_0_0_0, tg_yz_xxyyy_s_1_0_0, tg_yz_xxyyz_s_0_0_0, tg_yz_xxyyz_s_1_0_0, tg_yz_xxyzz_s_0_0_0, tg_yz_xxyzz_s_1_0_0, tg_yz_xxzzz_s_0_0_0, tg_yz_xxzzz_s_1_0_0, tg_yz_xyyyy_s_0_0_0, tg_yz_xyyyy_s_1_0_0, tg_yz_xyyyz_s_0_0_0, tg_yz_xyyyz_s_1_0_0, tg_yz_xyyzz_s_0_0_0, tg_yz_xyyzz_s_1_0_0, tg_yz_xyzzz_s_0_0_0, tg_yz_xyzzz_s_1_0_0, tg_yz_xzzzz_s_0_0_0, tg_yz_xzzzz_s_1_0_0, tg_yz_yyyyy_s_0_0_0, tg_yz_yyyyy_s_1_0_0, tg_yz_yyyyz_s_0_0_0, tg_yz_yyyyz_s_1_0_0, tg_yz_yyyzz_s_0_0_0, tg_yz_yyyzz_s_1_0_0, tg_yz_yyzzz_s_0_0_0, tg_yz_yyzzz_s_1_0_0, tg_yz_yzzzz_s_0_0_0, tg_yz_yzzzz_s_1_0_0, tg_yz_zzzzz_s_0_0_0, tg_yz_zzzzz_s_1_0_0, tg_yzz_xxxxx_s_0_0_0, tg_yzz_xxxxy_s_0_0_0, tg_yzz_xxxxz_s_0_0_0, tg_yzz_xxxyy_s_0_0_0, tg_yzz_xxxyz_s_0_0_0, tg_yzz_xxxzz_s_0_0_0, tg_yzz_xxyyy_s_0_0_0, tg_yzz_xxyyz_s_0_0_0, tg_yzz_xxyzz_s_0_0_0, tg_yzz_xxzzz_s_0_0_0, tg_yzz_xyyyy_s_0_0_0, tg_yzz_xyyyz_s_0_0_0, tg_yzz_xyyzz_s_0_0_0, tg_yzz_xyzzz_s_0_0_0, tg_yzz_xzzzz_s_0_0_0, tg_yzz_yyyyy_s_0_0_0, tg_yzz_yyyyz_s_0_0_0, tg_yzz_yyyzz_s_0_0_0, tg_yzz_yyzzz_s_0_0_0, tg_yzz_yzzzz_s_0_0_0, tg_yzz_zzzzz_s_0_0_0, tg_z_xxxxx_s_0_0_0, tg_z_xxxxx_s_1_0_0, tg_z_xxxxy_s_0_0_0, tg_z_xxxxy_s_1_0_0, tg_z_xxxxz_s_0_0_0, tg_z_xxxxz_s_1_0_0, tg_z_xxxyy_s_0_0_0, tg_z_xxxyy_s_1_0_0, tg_z_xxxyz_s_0_0_0, tg_z_xxxyz_s_1_0_0, tg_z_xxxzz_s_0_0_0, tg_z_xxxzz_s_1_0_0, tg_z_xxyyy_s_0_0_0, tg_z_xxyyy_s_1_0_0, tg_z_xxyyz_s_0_0_0, tg_z_xxyyz_s_1_0_0, tg_z_xxyzz_s_0_0_0, tg_z_xxyzz_s_1_0_0, tg_z_xxzzz_s_0_0_0, tg_z_xxzzz_s_1_0_0, tg_z_xyyyy_s_0_0_0, tg_z_xyyyy_s_1_0_0, tg_z_xyyyz_s_0_0_0, tg_z_xyyyz_s_1_0_0, tg_z_xyyzz_s_0_0_0, tg_z_xyyzz_s_1_0_0, tg_z_xyzzz_s_0_0_0, tg_z_xyzzz_s_1_0_0, tg_z_xzzzz_s_0_0_0, tg_z_xzzzz_s_1_0_0, tg_z_yyyyy_s_0_0_0, tg_z_yyyyy_s_1_0_0, tg_z_yyyyz_s_0_0_0, tg_z_yyyyz_s_1_0_0, tg_z_yyyzz_s_0_0_0, tg_z_yyyzz_s_1_0_0, tg_z_yyzzz_s_0_0_0, tg_z_yyzzz_s_1_0_0, tg_z_yzzzz_s_0_0_0, tg_z_yzzzz_s_1_0_0, tg_z_zzzzz_s_0_0_0, tg_z_zzzzz_s_1_0_0, tg_zz_xxxxx_s_0_0_0, tg_zz_xxxxx_s_1_0_0, tg_zz_xxxxy_s_0_0_0, tg_zz_xxxxy_s_1_0_0, tg_zz_xxxxz_s_0_0_0, tg_zz_xxxxz_s_1_0_0, tg_zz_xxxyy_s_0_0_0, tg_zz_xxxyy_s_1_0_0, tg_zz_xxxyz_s_0_0_0, tg_zz_xxxyz_s_1_0_0, tg_zz_xxxzz_s_0_0_0, tg_zz_xxxzz_s_1_0_0, tg_zz_xxyyy_s_0_0_0, tg_zz_xxyyy_s_1_0_0, tg_zz_xxyyz_s_0_0_0, tg_zz_xxyyz_s_1_0_0, tg_zz_xxyzz_s_0_0_0, tg_zz_xxyzz_s_1_0_0, tg_zz_xxzzz_s_0_0_0, tg_zz_xxzzz_s_1_0_0, tg_zz_xyyyy_s_0_0_0, tg_zz_xyyyy_s_1_0_0, tg_zz_xyyyz_s_0_0_0, tg_zz_xyyyz_s_1_0_0, tg_zz_xyyzz_s_0_0_0, tg_zz_xyyzz_s_1_0_0, tg_zz_xyzzz_s_0_0_0, tg_zz_xyzzz_s_1_0_0, tg_zz_xzzzz_s_0_0_0, tg_zz_xzzzz_s_1_0_0, tg_zz_yyyyy_s_0_0_0, tg_zz_yyyyy_s_1_0_0, tg_zz_yyyyz_s_0_0_0, tg_zz_yyyyz_s_1_0_0, tg_zz_yyyzz_s_0_0_0, tg_zz_yyyzz_s_1_0_0, tg_zz_yyzzz_s_0_0_0, tg_zz_yyzzz_s_1_0_0, tg_zz_yzzzz_s_0_0_0, tg_zz_yzzzz_s_1_0_0, tg_zz_zzzzz_s_0_0_0, tg_zz_zzzzz_s_1_0_0, tg_zzz_xxxxx_s_0_0_0, tg_zzz_xxxxy_s_0_0_0, tg_zzz_xxxxz_s_0_0_0, tg_zzz_xxxyy_s_0_0_0, tg_zzz_xxxyz_s_0_0_0, tg_zzz_xxxzz_s_0_0_0, tg_zzz_xxyyy_s_0_0_0, tg_zzz_xxyyz_s_0_0_0, tg_zzz_xxyzz_s_0_0_0, tg_zzz_xxzzz_s_0_0_0, tg_zzz_xyyyy_s_0_0_0, tg_zzz_xyyyz_s_0_0_0, tg_zzz_xyyzz_s_0_0_0, tg_zzz_xyzzz_s_0_0_0, tg_zzz_xzzzz_s_0_0_0, tg_zzz_yyyyy_s_0_0_0, tg_zzz_yyyyz_s_0_0_0, tg_zzz_yyyzz_s_0_0_0, tg_zzz_yyzzz_s_0_0_0, tg_zzz_yzzzz_s_0_0_0, tg_zzz_zzzzz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_xxx_xxxxx_s_0_0_0[i] = 2.0 * tg_x_xxxxx_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxxy_s_0_0_0[i] = 2.0 * tg_x_xxxxy_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxxz_s_0_0_0[i] = 2.0 * tg_x_xxxxz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxyy_s_0_0_0[i] = 2.0 * tg_x_xxxyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxyz_s_0_0_0[i] = 2.0 * tg_x_xxxyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxzz_s_0_0_0[i] = 2.0 * tg_x_xxxzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxyyy_s_0_0_0[i] = 2.0 * tg_x_xxyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxyyz_s_0_0_0[i] = 2.0 * tg_x_xxyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxyzz_s_0_0_0[i] = 2.0 * tg_x_xxyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxzzz_s_0_0_0[i] = 2.0 * tg_x_xxzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyyyy_s_0_0_0[i] = 2.0 * tg_x_xyyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyyyz_s_0_0_0[i] = 2.0 * tg_x_xyyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyyzz_s_0_0_0[i] = 2.0 * tg_x_xyyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyzzz_s_0_0_0[i] = 2.0 * tg_x_xyzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_xzzzz_s_0_0_0[i] = 2.0 * tg_x_xzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_xzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_xzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyyyy_s_0_0_0[i] = 2.0 * tg_x_yyyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_yyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyyyz_s_0_0_0[i] = 2.0 * tg_x_yyyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_yyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyyzz_s_0_0_0[i] = 2.0 * tg_x_yyyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_yyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyzzz_s_0_0_0[i] = 2.0 * tg_x_yyzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_yyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_yzzzz_s_0_0_0[i] = 2.0 * tg_x_yzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_yzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_yzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxx_zzzzz_s_0_0_0[i] = 2.0 * tg_x_zzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_x_zzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_xx_zzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xxy_xxxxx_s_0_0_0[i] = 2.0 * tg_xx_xxxxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxx_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxxy_s_0_0_0[i] = 2.0 * tg_xx_xxxxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxy_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxxz_s_0_0_0[i] = 2.0 * tg_xx_xxxxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxyy_s_0_0_0[i] = 2.0 * tg_xx_xxxyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxyy_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxyz_s_0_0_0[i] = 2.0 * tg_xx_xxxyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxyz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxzz_s_0_0_0[i] = 2.0 * tg_xx_xxxzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxyyy_s_0_0_0[i] = 2.0 * tg_xx_xxyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyyy_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxyyz_s_0_0_0[i] = 2.0 * tg_xx_xxyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyyz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxyzz_s_0_0_0[i] = 2.0 * tg_xx_xxyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxzzz_s_0_0_0[i] = 2.0 * tg_xx_xxzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxzzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyyyy_s_0_0_0[i] = 2.0 * tg_xx_xyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyyyz_s_0_0_0[i] = 2.0 * tg_xx_xyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyyzz_s_0_0_0[i] = 2.0 * tg_xx_xyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyzzz_s_0_0_0[i] = 2.0 * tg_xx_xyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_xzzzz_s_0_0_0[i] = 2.0 * tg_xx_xzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyyyy_s_0_0_0[i] = 2.0 * tg_xx_yyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyyyz_s_0_0_0[i] = 2.0 * tg_xx_yyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyyzz_s_0_0_0[i] = 2.0 * tg_xx_yyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyzzz_s_0_0_0[i] = 2.0 * tg_xx_yyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_yzzzz_s_0_0_0[i] = 2.0 * tg_xx_yzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxy_zzzzz_s_0_0_0[i] = 2.0 * tg_xx_zzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_xxz_xxxxx_s_0_0_0[i] = 2.0 * tg_xx_xxxxx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxx_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxxy_s_0_0_0[i] = 2.0 * tg_xx_xxxxy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxy_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxxz_s_0_0_0[i] = 2.0 * tg_xx_xxxxz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxxz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxyy_s_0_0_0[i] = 2.0 * tg_xx_xxxyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxyy_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxyz_s_0_0_0[i] = 2.0 * tg_xx_xxxyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxyz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxzz_s_0_0_0[i] = 2.0 * tg_xx_xxxzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxyyy_s_0_0_0[i] = 2.0 * tg_xx_xxyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyyy_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxyyz_s_0_0_0[i] = 2.0 * tg_xx_xxyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyyz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxyzz_s_0_0_0[i] = 2.0 * tg_xx_xxyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxzzz_s_0_0_0[i] = 2.0 * tg_xx_xxzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxzzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyyyy_s_0_0_0[i] = 2.0 * tg_xx_xyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyyyz_s_0_0_0[i] = 2.0 * tg_xx_xyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyyzz_s_0_0_0[i] = 2.0 * tg_xx_xyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyzzz_s_0_0_0[i] = 2.0 * tg_xx_xyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_xzzzz_s_0_0_0[i] = 2.0 * tg_xx_xzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyyyy_s_0_0_0[i] = 2.0 * tg_xx_yyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyyyz_s_0_0_0[i] = 2.0 * tg_xx_yyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyyzz_s_0_0_0[i] = 2.0 * tg_xx_yyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyzzz_s_0_0_0[i] = 2.0 * tg_xx_yyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_yzzzz_s_0_0_0[i] = 2.0 * tg_xx_yzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_xxz_zzzzz_s_0_0_0[i] = 2.0 * tg_xx_zzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_xyy_xxxxx_s_0_0_0[i] = 2.0 * tg_yy_xxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxxy_s_0_0_0[i] = 2.0 * tg_yy_xxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxxz_s_0_0_0[i] = 2.0 * tg_yy_xxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxyy_s_0_0_0[i] = 2.0 * tg_yy_xxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxyz_s_0_0_0[i] = 2.0 * tg_yy_xxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxzz_s_0_0_0[i] = 2.0 * tg_yy_xxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxyyy_s_0_0_0[i] = 2.0 * tg_yy_xxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxyyz_s_0_0_0[i] = 2.0 * tg_yy_xxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxyzz_s_0_0_0[i] = 2.0 * tg_yy_xxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxzzz_s_0_0_0[i] = 2.0 * tg_yy_xxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyyyy_s_0_0_0[i] = 2.0 * tg_yy_xyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyyyz_s_0_0_0[i] = 2.0 * tg_yy_xyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyyzz_s_0_0_0[i] = 2.0 * tg_yy_xyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyzzz_s_0_0_0[i] = 2.0 * tg_yy_xyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_xzzzz_s_0_0_0[i] = 2.0 * tg_yy_xzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyyyy_s_0_0_0[i] = 2.0 * tg_yy_yyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyyyz_s_0_0_0[i] = 2.0 * tg_yy_yyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyyzz_s_0_0_0[i] = 2.0 * tg_yy_yyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyzzz_s_0_0_0[i] = 2.0 * tg_yy_yyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_yzzzz_s_0_0_0[i] = 2.0 * tg_yy_yzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyy_zzzzz_s_0_0_0[i] = 2.0 * tg_yy_zzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxxx_s_0_0_0[i] = 2.0 * tg_yz_xxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxxy_s_0_0_0[i] = 2.0 * tg_yz_xxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxxz_s_0_0_0[i] = 2.0 * tg_yz_xxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxyy_s_0_0_0[i] = 2.0 * tg_yz_xxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxyz_s_0_0_0[i] = 2.0 * tg_yz_xxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxzz_s_0_0_0[i] = 2.0 * tg_yz_xxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxyyy_s_0_0_0[i] = 2.0 * tg_yz_xxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxyyz_s_0_0_0[i] = 2.0 * tg_yz_xxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxyzz_s_0_0_0[i] = 2.0 * tg_yz_xxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxzzz_s_0_0_0[i] = 2.0 * tg_yz_xxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xyyyy_s_0_0_0[i] = 2.0 * tg_yz_xyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xyyyz_s_0_0_0[i] = 2.0 * tg_yz_xyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xyyzz_s_0_0_0[i] = 2.0 * tg_yz_xyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xyzzz_s_0_0_0[i] = 2.0 * tg_yz_xyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_xzzzz_s_0_0_0[i] = 2.0 * tg_yz_xzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyyyy_s_0_0_0[i] = 2.0 * tg_yz_yyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyyyz_s_0_0_0[i] = 2.0 * tg_yz_yyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyyzz_s_0_0_0[i] = 2.0 * tg_yz_yyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyzzz_s_0_0_0[i] = 2.0 * tg_yz_yyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_yzzzz_s_0_0_0[i] = 2.0 * tg_yz_yzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xyz_zzzzz_s_0_0_0[i] = 2.0 * tg_yz_zzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_zzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxxx_s_0_0_0[i] = 2.0 * tg_zz_xxxxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxx_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxxy_s_0_0_0[i] = 2.0 * tg_zz_xxxxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxy_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxxz_s_0_0_0[i] = 2.0 * tg_zz_xxxxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxyy_s_0_0_0[i] = 2.0 * tg_zz_xxxyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxyy_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxyz_s_0_0_0[i] = 2.0 * tg_zz_xxxyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxyz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxzz_s_0_0_0[i] = 2.0 * tg_zz_xxxzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxyyy_s_0_0_0[i] = 2.0 * tg_zz_xxyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxyyz_s_0_0_0[i] = 2.0 * tg_zz_xxyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxyzz_s_0_0_0[i] = 2.0 * tg_zz_xxyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxzzz_s_0_0_0[i] = 2.0 * tg_zz_xxzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyyyy_s_0_0_0[i] = 2.0 * tg_zz_xyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyyyz_s_0_0_0[i] = 2.0 * tg_zz_xyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyyzz_s_0_0_0[i] = 2.0 * tg_zz_xyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyzzz_s_0_0_0[i] = 2.0 * tg_zz_xyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_xzzzz_s_0_0_0[i] = 2.0 * tg_zz_xzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyyyy_s_0_0_0[i] = 2.0 * tg_zz_yyyyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyyy_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyyyz_s_0_0_0[i] = 2.0 * tg_zz_yyyyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyyz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyyzz_s_0_0_0[i] = 2.0 * tg_zz_yyyzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyzzz_s_0_0_0[i] = 2.0 * tg_zz_yyzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_yzzzz_s_0_0_0[i] = 2.0 * tg_zz_yzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_xzz_zzzzz_s_0_0_0[i] = 2.0 * tg_zz_zzzzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzzzz_s_0_0_0[i] * a_x * faz_0;

        tg_yyy_xxxxx_s_0_0_0[i] = 2.0 * tg_y_xxxxx_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxxxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxx_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxxy_s_0_0_0[i] = 2.0 * tg_y_xxxxy_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxxxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxy_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxxz_s_0_0_0[i] = 2.0 * tg_y_xxxxz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxxxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxyy_s_0_0_0[i] = 2.0 * tg_y_xxxyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxxyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxyy_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxyz_s_0_0_0[i] = 2.0 * tg_y_xxxyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxxyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxyz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxzz_s_0_0_0[i] = 2.0 * tg_y_xxxzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxxzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxyyy_s_0_0_0[i] = 2.0 * tg_y_xxyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxyyz_s_0_0_0[i] = 2.0 * tg_y_xxyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxyzz_s_0_0_0[i] = 2.0 * tg_y_xxyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxzzz_s_0_0_0[i] = 2.0 * tg_y_xxzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xxzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyyyy_s_0_0_0[i] = 2.0 * tg_y_xyyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyyyz_s_0_0_0[i] = 2.0 * tg_y_xyyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyyzz_s_0_0_0[i] = 2.0 * tg_y_xyyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyzzz_s_0_0_0[i] = 2.0 * tg_y_xyzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_xzzzz_s_0_0_0[i] = 2.0 * tg_y_xzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_xzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_xzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyyyy_s_0_0_0[i] = 2.0 * tg_y_yyyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_yyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyyyz_s_0_0_0[i] = 2.0 * tg_y_yyyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_yyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyyzz_s_0_0_0[i] = 2.0 * tg_y_yyyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_yyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyzzz_s_0_0_0[i] = 2.0 * tg_y_yyzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_yyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_yzzzz_s_0_0_0[i] = 2.0 * tg_y_yzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_yzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_yzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyy_zzzzz_s_0_0_0[i] = 2.0 * tg_y_zzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_y_zzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_yy_zzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yyz_xxxxx_s_0_0_0[i] = 2.0 * tg_yy_xxxxx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxx_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxxy_s_0_0_0[i] = 2.0 * tg_yy_xxxxy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxy_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxxz_s_0_0_0[i] = 2.0 * tg_yy_xxxxz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxxz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxyy_s_0_0_0[i] = 2.0 * tg_yy_xxxyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxyy_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxyz_s_0_0_0[i] = 2.0 * tg_yy_xxxyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxyz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxzz_s_0_0_0[i] = 2.0 * tg_yy_xxxzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxyyy_s_0_0_0[i] = 2.0 * tg_yy_xxyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyyy_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxyyz_s_0_0_0[i] = 2.0 * tg_yy_xxyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyyz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxyzz_s_0_0_0[i] = 2.0 * tg_yy_xxyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxzzz_s_0_0_0[i] = 2.0 * tg_yy_xxzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxzzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyyyy_s_0_0_0[i] = 2.0 * tg_yy_xyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyyyz_s_0_0_0[i] = 2.0 * tg_yy_xyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyyzz_s_0_0_0[i] = 2.0 * tg_yy_xyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyzzz_s_0_0_0[i] = 2.0 * tg_yy_xyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_xzzzz_s_0_0_0[i] = 2.0 * tg_yy_xzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyyyy_s_0_0_0[i] = 2.0 * tg_yy_yyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyyyz_s_0_0_0[i] = 2.0 * tg_yy_yyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyyzz_s_0_0_0[i] = 2.0 * tg_yy_yyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyzzz_s_0_0_0[i] = 2.0 * tg_yy_yyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_yzzzz_s_0_0_0[i] = 2.0 * tg_yy_yzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_yyz_zzzzz_s_0_0_0[i] = 2.0 * tg_yy_zzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_yzz_xxxxx_s_0_0_0[i] = 2.0 * tg_zz_xxxxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxx_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxxy_s_0_0_0[i] = 2.0 * tg_zz_xxxxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxy_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxxz_s_0_0_0[i] = 2.0 * tg_zz_xxxxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxyy_s_0_0_0[i] = 2.0 * tg_zz_xxxyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxyy_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxyz_s_0_0_0[i] = 2.0 * tg_zz_xxxyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxyz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxzz_s_0_0_0[i] = 2.0 * tg_zz_xxxzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxyyy_s_0_0_0[i] = 2.0 * tg_zz_xxyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxyyz_s_0_0_0[i] = 2.0 * tg_zz_xxyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxyzz_s_0_0_0[i] = 2.0 * tg_zz_xxyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxzzz_s_0_0_0[i] = 2.0 * tg_zz_xxzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyyyy_s_0_0_0[i] = 2.0 * tg_zz_xyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyyyz_s_0_0_0[i] = 2.0 * tg_zz_xyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyyzz_s_0_0_0[i] = 2.0 * tg_zz_xyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyzzz_s_0_0_0[i] = 2.0 * tg_zz_xyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_xzzzz_s_0_0_0[i] = 2.0 * tg_zz_xzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyyyy_s_0_0_0[i] = 2.0 * tg_zz_yyyyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyyy_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyyyz_s_0_0_0[i] = 2.0 * tg_zz_yyyyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyyz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyyzz_s_0_0_0[i] = 2.0 * tg_zz_yyyzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyzzz_s_0_0_0[i] = 2.0 * tg_zz_yyzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_yzzzz_s_0_0_0[i] = 2.0 * tg_zz_yzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_yzz_zzzzz_s_0_0_0[i] = 2.0 * tg_zz_zzzzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzzzz_s_0_0_0[i] * a_y * faz_0;

        tg_zzz_xxxxx_s_0_0_0[i] = 2.0 * tg_z_xxxxx_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxxxx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxx_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxxy_s_0_0_0[i] = 2.0 * tg_z_xxxxy_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxxxy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxy_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxxz_s_0_0_0[i] = 2.0 * tg_z_xxxxz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxxxz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxxz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxyy_s_0_0_0[i] = 2.0 * tg_z_xxxyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxxyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxyy_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxyz_s_0_0_0[i] = 2.0 * tg_z_xxxyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxxyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxyz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxzz_s_0_0_0[i] = 2.0 * tg_z_xxxzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxxzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxyyy_s_0_0_0[i] = 2.0 * tg_z_xxyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxyyz_s_0_0_0[i] = 2.0 * tg_z_xxyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxyzz_s_0_0_0[i] = 2.0 * tg_z_xxyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxzzz_s_0_0_0[i] = 2.0 * tg_z_xxzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xxzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyyyy_s_0_0_0[i] = 2.0 * tg_z_xyyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyyyz_s_0_0_0[i] = 2.0 * tg_z_xyyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyyzz_s_0_0_0[i] = 2.0 * tg_z_xyyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyzzz_s_0_0_0[i] = 2.0 * tg_z_xyzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_xzzzz_s_0_0_0[i] = 2.0 * tg_z_xzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_xzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_xzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyyyy_s_0_0_0[i] = 2.0 * tg_z_yyyyy_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyyyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_yyyyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyyy_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyyyz_s_0_0_0[i] = 2.0 * tg_z_yyyyz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyyyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_yyyyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyyz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyyzz_s_0_0_0[i] = 2.0 * tg_z_yyyzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyyzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_yyyzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyzzz_s_0_0_0[i] = 2.0 * tg_z_yyzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_yyzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_yzzzz_s_0_0_0[i] = 2.0 * tg_z_yzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_yzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_yzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzzzz_s_0_0_0[i] * a_z * faz_0;

        tg_zzz_zzzzz_s_0_0_0[i] = 2.0 * tg_z_zzzzz_s_0_0_0[i] * fzi_0 + 2.0 * tg_z_zzzzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_zz_zzzzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzzzz_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : PH

        auto tg_x_xxxxx_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1);

        auto tg_x_xxxxy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 1);

        auto tg_x_xxxxz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 2);

        auto tg_x_xxxyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 3);

        auto tg_x_xxxyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 4);

        auto tg_x_xxxzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 5);

        auto tg_x_xxyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 6);

        auto tg_x_xxyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 7);

        auto tg_x_xxyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 8);

        auto tg_x_xxzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 9);

        auto tg_x_xyyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 10);

        auto tg_x_xyyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 11);

        auto tg_x_xyyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 12);

        auto tg_x_xyzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 13);

        auto tg_x_xzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 14);

        auto tg_x_yyyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 15);

        auto tg_x_yyyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 16);

        auto tg_x_yyyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 17);

        auto tg_x_yyzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 18);

        auto tg_x_yzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 19);

        auto tg_x_zzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 20);

        auto tg_y_xxxxx_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 21);

        auto tg_y_xxxxy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 22);

        auto tg_y_xxxxz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 23);

        auto tg_y_xxxyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 24);

        auto tg_y_xxxyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 25);

        auto tg_y_xxxzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 26);

        auto tg_y_xxyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 27);

        auto tg_y_xxyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 28);

        auto tg_y_xxyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 29);

        auto tg_y_xxzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 30);

        auto tg_y_xyyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 31);

        auto tg_y_xyyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 32);

        auto tg_y_xyyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 33);

        auto tg_y_xyzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 34);

        auto tg_y_xzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 35);

        auto tg_y_yyyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 36);

        auto tg_y_yyyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 37);

        auto tg_y_yyyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 38);

        auto tg_y_yyzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 39);

        auto tg_y_yzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 40);

        auto tg_y_zzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 41);

        auto tg_z_xxxxx_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 42);

        auto tg_z_xxxxy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 43);

        auto tg_z_xxxxz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 44);

        auto tg_z_xxxyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 45);

        auto tg_z_xxxyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 46);

        auto tg_z_xxxzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 47);

        auto tg_z_xxyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 48);

        auto tg_z_xxyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 49);

        auto tg_z_xxyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 50);

        auto tg_z_xxzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 51);

        auto tg_z_xyyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 52);

        auto tg_z_xyyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 53);

        auto tg_z_xyyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 54);

        auto tg_z_xyzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 55);

        auto tg_z_xzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 56);

        auto tg_z_yyyyy_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 57);

        auto tg_z_yyyyz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 58);

        auto tg_z_yyyzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 59);

        auto tg_z_yyzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 60);

        auto tg_z_yzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 61);

        auto tg_z_zzzzz_s_0_0_1 = pbuffer.data(idx_ph_s_0_0_1 + 62);

        // Set up components of auxiliary buffer : DH

        auto tg_xx_xxxxx_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1);

        auto tg_xx_xxxxy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 1);

        auto tg_xx_xxxxz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 2);

        auto tg_xx_xxxyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 3);

        auto tg_xx_xxxyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 4);

        auto tg_xx_xxxzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 5);

        auto tg_xx_xxyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 6);

        auto tg_xx_xxyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 7);

        auto tg_xx_xxyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 8);

        auto tg_xx_xxzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 9);

        auto tg_xx_xyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 10);

        auto tg_xx_xyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 11);

        auto tg_xx_xyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 12);

        auto tg_xx_xyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 13);

        auto tg_xx_xzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 14);

        auto tg_xx_yyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 15);

        auto tg_xx_yyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 16);

        auto tg_xx_yyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 17);

        auto tg_xx_yyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 18);

        auto tg_xx_yzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 19);

        auto tg_xx_zzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 20);

        auto tg_xy_xxxxx_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 21);

        auto tg_xy_xxxxy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 22);

        auto tg_xy_xxxxz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 23);

        auto tg_xy_xxxyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 24);

        auto tg_xy_xxxyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 25);

        auto tg_xy_xxxzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 26);

        auto tg_xy_xxyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 27);

        auto tg_xy_xxyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 28);

        auto tg_xy_xxyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 29);

        auto tg_xy_xxzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 30);

        auto tg_xy_xyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 31);

        auto tg_xy_xyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 32);

        auto tg_xy_xyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 33);

        auto tg_xy_xyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 34);

        auto tg_xy_xzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 35);

        auto tg_xy_yyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 36);

        auto tg_xy_yyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 37);

        auto tg_xy_yyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 38);

        auto tg_xy_yyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 39);

        auto tg_xy_yzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 40);

        auto tg_xy_zzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 41);

        auto tg_xz_xxxxx_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 42);

        auto tg_xz_xxxxy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 43);

        auto tg_xz_xxxxz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 44);

        auto tg_xz_xxxyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 45);

        auto tg_xz_xxxyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 46);

        auto tg_xz_xxxzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 47);

        auto tg_xz_xxyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 48);

        auto tg_xz_xxyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 49);

        auto tg_xz_xxyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 50);

        auto tg_xz_xxzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 51);

        auto tg_xz_xyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 52);

        auto tg_xz_xyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 53);

        auto tg_xz_xyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 54);

        auto tg_xz_xyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 55);

        auto tg_xz_xzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 56);

        auto tg_xz_yyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 57);

        auto tg_xz_yyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 58);

        auto tg_xz_yyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 59);

        auto tg_xz_yyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 60);

        auto tg_xz_yzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 61);

        auto tg_xz_zzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 62);

        auto tg_yy_xxxxx_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 63);

        auto tg_yy_xxxxy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 64);

        auto tg_yy_xxxxz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 65);

        auto tg_yy_xxxyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 66);

        auto tg_yy_xxxyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 67);

        auto tg_yy_xxxzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 68);

        auto tg_yy_xxyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 69);

        auto tg_yy_xxyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 70);

        auto tg_yy_xxyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 71);

        auto tg_yy_xxzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 72);

        auto tg_yy_xyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 73);

        auto tg_yy_xyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 74);

        auto tg_yy_xyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 75);

        auto tg_yy_xyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 76);

        auto tg_yy_xzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 77);

        auto tg_yy_yyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 78);

        auto tg_yy_yyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 79);

        auto tg_yy_yyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 80);

        auto tg_yy_yyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 81);

        auto tg_yy_yzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 82);

        auto tg_yy_zzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 83);

        auto tg_yz_xxxxx_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 84);

        auto tg_yz_xxxxy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 85);

        auto tg_yz_xxxxz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 86);

        auto tg_yz_xxxyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 87);

        auto tg_yz_xxxyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 88);

        auto tg_yz_xxxzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 89);

        auto tg_yz_xxyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 90);

        auto tg_yz_xxyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 91);

        auto tg_yz_xxyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 92);

        auto tg_yz_xxzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 93);

        auto tg_yz_xyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 94);

        auto tg_yz_xyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 95);

        auto tg_yz_xyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 96);

        auto tg_yz_xyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 97);

        auto tg_yz_xzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 98);

        auto tg_yz_yyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 99);

        auto tg_yz_yyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 100);

        auto tg_yz_yyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 101);

        auto tg_yz_yyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 102);

        auto tg_yz_yzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 103);

        auto tg_yz_zzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 104);

        auto tg_zz_xxxxx_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 105);

        auto tg_zz_xxxxy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 106);

        auto tg_zz_xxxxz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 107);

        auto tg_zz_xxxyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 108);

        auto tg_zz_xxxyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 109);

        auto tg_zz_xxxzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 110);

        auto tg_zz_xxyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 111);

        auto tg_zz_xxyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 112);

        auto tg_zz_xxyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 113);

        auto tg_zz_xxzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 114);

        auto tg_zz_xyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 115);

        auto tg_zz_xyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 116);

        auto tg_zz_xyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 117);

        auto tg_zz_xyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 118);

        auto tg_zz_xzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 119);

        auto tg_zz_yyyyy_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 120);

        auto tg_zz_yyyyz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 121);

        auto tg_zz_yyyzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 122);

        auto tg_zz_yyzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 123);

        auto tg_zz_yzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 124);

        auto tg_zz_zzzzz_s_0_0_1 = pbuffer.data(idx_dh_s_0_0_1 + 125);

        #pragma omp simd aligned(b_exps, tg_x_xxxxx_s_0_0_1, tg_x_xxxxy_s_0_0_1, tg_x_xxxxz_s_0_0_1, tg_x_xxxyy_s_0_0_1, tg_x_xxxyz_s_0_0_1, tg_x_xxxzz_s_0_0_1, tg_x_xxyyy_s_0_0_1, tg_x_xxyyz_s_0_0_1, tg_x_xxyzz_s_0_0_1, tg_x_xxzzz_s_0_0_1, tg_x_xyyyy_s_0_0_1, tg_x_xyyyz_s_0_0_1, tg_x_xyyzz_s_0_0_1, tg_x_xyzzz_s_0_0_1, tg_x_xzzzz_s_0_0_1, tg_x_yyyyy_s_0_0_1, tg_x_yyyyz_s_0_0_1, tg_x_yyyzz_s_0_0_1, tg_x_yyzzz_s_0_0_1, tg_x_yzzzz_s_0_0_1, tg_x_zzzzz_s_0_0_1, tg_xx_xxxxx_s_0_0_1, tg_xx_xxxxy_s_0_0_1, tg_xx_xxxxz_s_0_0_1, tg_xx_xxxyy_s_0_0_1, tg_xx_xxxyz_s_0_0_1, tg_xx_xxxzz_s_0_0_1, tg_xx_xxyyy_s_0_0_1, tg_xx_xxyyz_s_0_0_1, tg_xx_xxyzz_s_0_0_1, tg_xx_xxzzz_s_0_0_1, tg_xx_xyyyy_s_0_0_1, tg_xx_xyyyz_s_0_0_1, tg_xx_xyyzz_s_0_0_1, tg_xx_xyzzz_s_0_0_1, tg_xx_xzzzz_s_0_0_1, tg_xx_yyyyy_s_0_0_1, tg_xx_yyyyz_s_0_0_1, tg_xx_yyyzz_s_0_0_1, tg_xx_yyzzz_s_0_0_1, tg_xx_yzzzz_s_0_0_1, tg_xx_zzzzz_s_0_0_1, tg_xxx_xxxxx_s_0_0_0, tg_xxx_xxxxy_s_0_0_0, tg_xxx_xxxxz_s_0_0_0, tg_xxx_xxxyy_s_0_0_0, tg_xxx_xxxyz_s_0_0_0, tg_xxx_xxxzz_s_0_0_0, tg_xxx_xxyyy_s_0_0_0, tg_xxx_xxyyz_s_0_0_0, tg_xxx_xxyzz_s_0_0_0, tg_xxx_xxzzz_s_0_0_0, tg_xxx_xyyyy_s_0_0_0, tg_xxx_xyyyz_s_0_0_0, tg_xxx_xyyzz_s_0_0_0, tg_xxx_xyzzz_s_0_0_0, tg_xxx_xzzzz_s_0_0_0, tg_xxx_yyyyy_s_0_0_0, tg_xxx_yyyyz_s_0_0_0, tg_xxx_yyyzz_s_0_0_0, tg_xxx_yyzzz_s_0_0_0, tg_xxx_yzzzz_s_0_0_0, tg_xxx_zzzzz_s_0_0_0, tg_xxy_xxxxx_s_0_0_0, tg_xxy_xxxxy_s_0_0_0, tg_xxy_xxxxz_s_0_0_0, tg_xxy_xxxyy_s_0_0_0, tg_xxy_xxxyz_s_0_0_0, tg_xxy_xxxzz_s_0_0_0, tg_xxy_xxyyy_s_0_0_0, tg_xxy_xxyyz_s_0_0_0, tg_xxy_xxyzz_s_0_0_0, tg_xxy_xxzzz_s_0_0_0, tg_xxy_xyyyy_s_0_0_0, tg_xxy_xyyyz_s_0_0_0, tg_xxy_xyyzz_s_0_0_0, tg_xxy_xyzzz_s_0_0_0, tg_xxy_xzzzz_s_0_0_0, tg_xxy_yyyyy_s_0_0_0, tg_xxy_yyyyz_s_0_0_0, tg_xxy_yyyzz_s_0_0_0, tg_xxy_yyzzz_s_0_0_0, tg_xxy_yzzzz_s_0_0_0, tg_xxy_zzzzz_s_0_0_0, tg_xxz_xxxxx_s_0_0_0, tg_xxz_xxxxy_s_0_0_0, tg_xxz_xxxxz_s_0_0_0, tg_xxz_xxxyy_s_0_0_0, tg_xxz_xxxyz_s_0_0_0, tg_xxz_xxxzz_s_0_0_0, tg_xxz_xxyyy_s_0_0_0, tg_xxz_xxyyz_s_0_0_0, tg_xxz_xxyzz_s_0_0_0, tg_xxz_xxzzz_s_0_0_0, tg_xxz_xyyyy_s_0_0_0, tg_xxz_xyyyz_s_0_0_0, tg_xxz_xyyzz_s_0_0_0, tg_xxz_xyzzz_s_0_0_0, tg_xxz_xzzzz_s_0_0_0, tg_xxz_yyyyy_s_0_0_0, tg_xxz_yyyyz_s_0_0_0, tg_xxz_yyyzz_s_0_0_0, tg_xxz_yyzzz_s_0_0_0, tg_xxz_yzzzz_s_0_0_0, tg_xxz_zzzzz_s_0_0_0, tg_xyy_xxxxx_s_0_0_0, tg_xyy_xxxxy_s_0_0_0, tg_xyy_xxxxz_s_0_0_0, tg_xyy_xxxyy_s_0_0_0, tg_xyy_xxxyz_s_0_0_0, tg_xyy_xxxzz_s_0_0_0, tg_xyy_xxyyy_s_0_0_0, tg_xyy_xxyyz_s_0_0_0, tg_xyy_xxyzz_s_0_0_0, tg_xyy_xxzzz_s_0_0_0, tg_xyy_xyyyy_s_0_0_0, tg_xyy_xyyyz_s_0_0_0, tg_xyy_xyyzz_s_0_0_0, tg_xyy_xyzzz_s_0_0_0, tg_xyy_xzzzz_s_0_0_0, tg_xyy_yyyyy_s_0_0_0, tg_xyy_yyyyz_s_0_0_0, tg_xyy_yyyzz_s_0_0_0, tg_xyy_yyzzz_s_0_0_0, tg_xyy_yzzzz_s_0_0_0, tg_xyy_zzzzz_s_0_0_0, tg_xyz_xxxxx_s_0_0_0, tg_xyz_xxxxy_s_0_0_0, tg_xyz_xxxxz_s_0_0_0, tg_xyz_xxxyy_s_0_0_0, tg_xyz_xxxyz_s_0_0_0, tg_xyz_xxxzz_s_0_0_0, tg_xyz_xxyyy_s_0_0_0, tg_xyz_xxyyz_s_0_0_0, tg_xyz_xxyzz_s_0_0_0, tg_xyz_xxzzz_s_0_0_0, tg_xyz_xyyyy_s_0_0_0, tg_xyz_xyyyz_s_0_0_0, tg_xyz_xyyzz_s_0_0_0, tg_xyz_xyzzz_s_0_0_0, tg_xyz_xzzzz_s_0_0_0, tg_xyz_yyyyy_s_0_0_0, tg_xyz_yyyyz_s_0_0_0, tg_xyz_yyyzz_s_0_0_0, tg_xyz_yyzzz_s_0_0_0, tg_xyz_yzzzz_s_0_0_0, tg_xyz_zzzzz_s_0_0_0, tg_xzz_xxxxx_s_0_0_0, tg_xzz_xxxxy_s_0_0_0, tg_xzz_xxxxz_s_0_0_0, tg_xzz_xxxyy_s_0_0_0, tg_xzz_xxxyz_s_0_0_0, tg_xzz_xxxzz_s_0_0_0, tg_xzz_xxyyy_s_0_0_0, tg_xzz_xxyyz_s_0_0_0, tg_xzz_xxyzz_s_0_0_0, tg_xzz_xxzzz_s_0_0_0, tg_xzz_xyyyy_s_0_0_0, tg_xzz_xyyyz_s_0_0_0, tg_xzz_xyyzz_s_0_0_0, tg_xzz_xyzzz_s_0_0_0, tg_xzz_xzzzz_s_0_0_0, tg_xzz_yyyyy_s_0_0_0, tg_xzz_yyyyz_s_0_0_0, tg_xzz_yyyzz_s_0_0_0, tg_xzz_yyzzz_s_0_0_0, tg_xzz_yzzzz_s_0_0_0, tg_xzz_zzzzz_s_0_0_0, tg_y_xxxxx_s_0_0_1, tg_y_xxxxy_s_0_0_1, tg_y_xxxxz_s_0_0_1, tg_y_xxxyy_s_0_0_1, tg_y_xxxyz_s_0_0_1, tg_y_xxxzz_s_0_0_1, tg_y_xxyyy_s_0_0_1, tg_y_xxyyz_s_0_0_1, tg_y_xxyzz_s_0_0_1, tg_y_xxzzz_s_0_0_1, tg_y_xyyyy_s_0_0_1, tg_y_xyyyz_s_0_0_1, tg_y_xyyzz_s_0_0_1, tg_y_xyzzz_s_0_0_1, tg_y_xzzzz_s_0_0_1, tg_y_yyyyy_s_0_0_1, tg_y_yyyyz_s_0_0_1, tg_y_yyyzz_s_0_0_1, tg_y_yyzzz_s_0_0_1, tg_y_yzzzz_s_0_0_1, tg_y_zzzzz_s_0_0_1, tg_yy_xxxxx_s_0_0_1, tg_yy_xxxxy_s_0_0_1, tg_yy_xxxxz_s_0_0_1, tg_yy_xxxyy_s_0_0_1, tg_yy_xxxyz_s_0_0_1, tg_yy_xxxzz_s_0_0_1, tg_yy_xxyyy_s_0_0_1, tg_yy_xxyyz_s_0_0_1, tg_yy_xxyzz_s_0_0_1, tg_yy_xxzzz_s_0_0_1, tg_yy_xyyyy_s_0_0_1, tg_yy_xyyyz_s_0_0_1, tg_yy_xyyzz_s_0_0_1, tg_yy_xyzzz_s_0_0_1, tg_yy_xzzzz_s_0_0_1, tg_yy_yyyyy_s_0_0_1, tg_yy_yyyyz_s_0_0_1, tg_yy_yyyzz_s_0_0_1, tg_yy_yyzzz_s_0_0_1, tg_yy_yzzzz_s_0_0_1, tg_yy_zzzzz_s_0_0_1, tg_yyy_xxxxx_s_0_0_0, tg_yyy_xxxxy_s_0_0_0, tg_yyy_xxxxz_s_0_0_0, tg_yyy_xxxyy_s_0_0_0, tg_yyy_xxxyz_s_0_0_0, tg_yyy_xxxzz_s_0_0_0, tg_yyy_xxyyy_s_0_0_0, tg_yyy_xxyyz_s_0_0_0, tg_yyy_xxyzz_s_0_0_0, tg_yyy_xxzzz_s_0_0_0, tg_yyy_xyyyy_s_0_0_0, tg_yyy_xyyyz_s_0_0_0, tg_yyy_xyyzz_s_0_0_0, tg_yyy_xyzzz_s_0_0_0, tg_yyy_xzzzz_s_0_0_0, tg_yyy_yyyyy_s_0_0_0, tg_yyy_yyyyz_s_0_0_0, tg_yyy_yyyzz_s_0_0_0, tg_yyy_yyzzz_s_0_0_0, tg_yyy_yzzzz_s_0_0_0, tg_yyy_zzzzz_s_0_0_0, tg_yyz_xxxxx_s_0_0_0, tg_yyz_xxxxy_s_0_0_0, tg_yyz_xxxxz_s_0_0_0, tg_yyz_xxxyy_s_0_0_0, tg_yyz_xxxyz_s_0_0_0, tg_yyz_xxxzz_s_0_0_0, tg_yyz_xxyyy_s_0_0_0, tg_yyz_xxyyz_s_0_0_0, tg_yyz_xxyzz_s_0_0_0, tg_yyz_xxzzz_s_0_0_0, tg_yyz_xyyyy_s_0_0_0, tg_yyz_xyyyz_s_0_0_0, tg_yyz_xyyzz_s_0_0_0, tg_yyz_xyzzz_s_0_0_0, tg_yyz_xzzzz_s_0_0_0, tg_yyz_yyyyy_s_0_0_0, tg_yyz_yyyyz_s_0_0_0, tg_yyz_yyyzz_s_0_0_0, tg_yyz_yyzzz_s_0_0_0, tg_yyz_yzzzz_s_0_0_0, tg_yyz_zzzzz_s_0_0_0, tg_yz_xxxxx_s_0_0_1, tg_yz_xxxxy_s_0_0_1, tg_yz_xxxxz_s_0_0_1, tg_yz_xxxyy_s_0_0_1, tg_yz_xxxyz_s_0_0_1, tg_yz_xxxzz_s_0_0_1, tg_yz_xxyyy_s_0_0_1, tg_yz_xxyyz_s_0_0_1, tg_yz_xxyzz_s_0_0_1, tg_yz_xxzzz_s_0_0_1, tg_yz_xyyyy_s_0_0_1, tg_yz_xyyyz_s_0_0_1, tg_yz_xyyzz_s_0_0_1, tg_yz_xyzzz_s_0_0_1, tg_yz_xzzzz_s_0_0_1, tg_yz_yyyyy_s_0_0_1, tg_yz_yyyyz_s_0_0_1, tg_yz_yyyzz_s_0_0_1, tg_yz_yyzzz_s_0_0_1, tg_yz_yzzzz_s_0_0_1, tg_yz_zzzzz_s_0_0_1, tg_yzz_xxxxx_s_0_0_0, tg_yzz_xxxxy_s_0_0_0, tg_yzz_xxxxz_s_0_0_0, tg_yzz_xxxyy_s_0_0_0, tg_yzz_xxxyz_s_0_0_0, tg_yzz_xxxzz_s_0_0_0, tg_yzz_xxyyy_s_0_0_0, tg_yzz_xxyyz_s_0_0_0, tg_yzz_xxyzz_s_0_0_0, tg_yzz_xxzzz_s_0_0_0, tg_yzz_xyyyy_s_0_0_0, tg_yzz_xyyyz_s_0_0_0, tg_yzz_xyyzz_s_0_0_0, tg_yzz_xyzzz_s_0_0_0, tg_yzz_xzzzz_s_0_0_0, tg_yzz_yyyyy_s_0_0_0, tg_yzz_yyyyz_s_0_0_0, tg_yzz_yyyzz_s_0_0_0, tg_yzz_yyzzz_s_0_0_0, tg_yzz_yzzzz_s_0_0_0, tg_yzz_zzzzz_s_0_0_0, tg_z_xxxxx_s_0_0_1, tg_z_xxxxy_s_0_0_1, tg_z_xxxxz_s_0_0_1, tg_z_xxxyy_s_0_0_1, tg_z_xxxyz_s_0_0_1, tg_z_xxxzz_s_0_0_1, tg_z_xxyyy_s_0_0_1, tg_z_xxyyz_s_0_0_1, tg_z_xxyzz_s_0_0_1, tg_z_xxzzz_s_0_0_1, tg_z_xyyyy_s_0_0_1, tg_z_xyyyz_s_0_0_1, tg_z_xyyzz_s_0_0_1, tg_z_xyzzz_s_0_0_1, tg_z_xzzzz_s_0_0_1, tg_z_yyyyy_s_0_0_1, tg_z_yyyyz_s_0_0_1, tg_z_yyyzz_s_0_0_1, tg_z_yyzzz_s_0_0_1, tg_z_yzzzz_s_0_0_1, tg_z_zzzzz_s_0_0_1, tg_zz_xxxxx_s_0_0_1, tg_zz_xxxxy_s_0_0_1, tg_zz_xxxxz_s_0_0_1, tg_zz_xxxyy_s_0_0_1, tg_zz_xxxyz_s_0_0_1, tg_zz_xxxzz_s_0_0_1, tg_zz_xxyyy_s_0_0_1, tg_zz_xxyyz_s_0_0_1, tg_zz_xxyzz_s_0_0_1, tg_zz_xxzzz_s_0_0_1, tg_zz_xyyyy_s_0_0_1, tg_zz_xyyyz_s_0_0_1, tg_zz_xyyzz_s_0_0_1, tg_zz_xyzzz_s_0_0_1, tg_zz_xzzzz_s_0_0_1, tg_zz_yyyyy_s_0_0_1, tg_zz_yyyyz_s_0_0_1, tg_zz_yyyzz_s_0_0_1, tg_zz_yyzzz_s_0_0_1, tg_zz_yzzzz_s_0_0_1, tg_zz_zzzzz_s_0_0_1, tg_zzz_xxxxx_s_0_0_0, tg_zzz_xxxxy_s_0_0_0, tg_zzz_xxxxz_s_0_0_0, tg_zzz_xxxyy_s_0_0_0, tg_zzz_xxxyz_s_0_0_0, tg_zzz_xxxzz_s_0_0_0, tg_zzz_xxyyy_s_0_0_0, tg_zzz_xxyyz_s_0_0_0, tg_zzz_xxyzz_s_0_0_0, tg_zzz_xxzzz_s_0_0_0, tg_zzz_xyyyy_s_0_0_0, tg_zzz_xyyyz_s_0_0_0, tg_zzz_xyyzz_s_0_0_0, tg_zzz_xyzzz_s_0_0_0, tg_zzz_xzzzz_s_0_0_0, tg_zzz_yyyyy_s_0_0_0, tg_zzz_yyyyz_s_0_0_0, tg_zzz_yyyzz_s_0_0_0, tg_zzz_yyzzz_s_0_0_0, tg_zzz_yzzzz_s_0_0_0, tg_zzz_zzzzz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxx_xxxxx_s_0_0_0[i] += tg_x_xxxxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxxy_s_0_0_0[i] += tg_x_xxxxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxxz_s_0_0_0[i] += tg_x_xxxxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxyy_s_0_0_0[i] += tg_x_xxxyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxyz_s_0_0_0[i] += tg_x_xxxyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxzz_s_0_0_0[i] += tg_x_xxxzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxyyy_s_0_0_0[i] += tg_x_xxyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxyyz_s_0_0_0[i] += tg_x_xxyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxyzz_s_0_0_0[i] += tg_x_xxyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxzzz_s_0_0_0[i] += tg_x_xxzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyyyy_s_0_0_0[i] += tg_x_xyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyyyz_s_0_0_0[i] += tg_x_xyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyyzz_s_0_0_0[i] += tg_x_xyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyzzz_s_0_0_0[i] += tg_x_xyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xzzzz_s_0_0_0[i] += tg_x_xzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyyyy_s_0_0_0[i] += tg_x_yyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyyyz_s_0_0_0[i] += tg_x_yyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyyzz_s_0_0_0[i] += tg_x_yyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyzzz_s_0_0_0[i] += tg_x_yyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yzzzz_s_0_0_0[i] += tg_x_yzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_zzzzz_s_0_0_0[i] += tg_x_zzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_zzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxy_xxxxx_s_0_0_0[i] += tg_xx_xxxxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxxy_s_0_0_0[i] += tg_xx_xxxxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxxz_s_0_0_0[i] += tg_xx_xxxxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxyy_s_0_0_0[i] += tg_xx_xxxyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxyz_s_0_0_0[i] += tg_xx_xxxyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxzz_s_0_0_0[i] += tg_xx_xxxzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxyyy_s_0_0_0[i] += tg_xx_xxyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxyyz_s_0_0_0[i] += tg_xx_xxyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxyzz_s_0_0_0[i] += tg_xx_xxyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxzzz_s_0_0_0[i] += tg_xx_xxzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyyyy_s_0_0_0[i] += tg_xx_xyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyyyz_s_0_0_0[i] += tg_xx_xyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyyzz_s_0_0_0[i] += tg_xx_xyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyzzz_s_0_0_0[i] += tg_xx_xyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xzzzz_s_0_0_0[i] += tg_xx_xzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyyyy_s_0_0_0[i] += tg_xx_yyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyyyz_s_0_0_0[i] += tg_xx_yyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyyzz_s_0_0_0[i] += tg_xx_yyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyzzz_s_0_0_0[i] += tg_xx_yyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yzzzz_s_0_0_0[i] += tg_xx_yzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_zzzzz_s_0_0_0[i] += tg_xx_zzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxz_xxxxx_s_0_0_0[i] += tg_xx_xxxxx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxxy_s_0_0_0[i] += tg_xx_xxxxy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxxz_s_0_0_0[i] += tg_xx_xxxxz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxyy_s_0_0_0[i] += tg_xx_xxxyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxyz_s_0_0_0[i] += tg_xx_xxxyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxzz_s_0_0_0[i] += tg_xx_xxxzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxyyy_s_0_0_0[i] += tg_xx_xxyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxyyz_s_0_0_0[i] += tg_xx_xxyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxyzz_s_0_0_0[i] += tg_xx_xxyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxzzz_s_0_0_0[i] += tg_xx_xxzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyyyy_s_0_0_0[i] += tg_xx_xyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyyyz_s_0_0_0[i] += tg_xx_xyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyyzz_s_0_0_0[i] += tg_xx_xyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyzzz_s_0_0_0[i] += tg_xx_xyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xzzzz_s_0_0_0[i] += tg_xx_xzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyyyy_s_0_0_0[i] += tg_xx_yyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyyyz_s_0_0_0[i] += tg_xx_yyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyyzz_s_0_0_0[i] += tg_xx_yyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyzzz_s_0_0_0[i] += tg_xx_yyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yzzzz_s_0_0_0[i] += tg_xx_yzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_zzzzz_s_0_0_0[i] += tg_xx_zzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xyy_xxxxx_s_0_0_0[i] += tg_yy_xxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxxy_s_0_0_0[i] += tg_yy_xxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxxz_s_0_0_0[i] += tg_yy_xxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxyy_s_0_0_0[i] += tg_yy_xxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxyz_s_0_0_0[i] += tg_yy_xxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxzz_s_0_0_0[i] += tg_yy_xxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxyyy_s_0_0_0[i] += tg_yy_xxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxyyz_s_0_0_0[i] += tg_yy_xxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxyzz_s_0_0_0[i] += tg_yy_xxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxzzz_s_0_0_0[i] += tg_yy_xxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyyyy_s_0_0_0[i] += tg_yy_xyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyyyz_s_0_0_0[i] += tg_yy_xyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyyzz_s_0_0_0[i] += tg_yy_xyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyzzz_s_0_0_0[i] += tg_yy_xyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xzzzz_s_0_0_0[i] += tg_yy_xzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyyyy_s_0_0_0[i] += tg_yy_yyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyyyz_s_0_0_0[i] += tg_yy_yyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyyzz_s_0_0_0[i] += tg_yy_yyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyzzz_s_0_0_0[i] += tg_yy_yyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yzzzz_s_0_0_0[i] += tg_yy_yzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_zzzzz_s_0_0_0[i] += tg_yy_zzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxxx_s_0_0_0[i] += tg_yz_xxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxxy_s_0_0_0[i] += tg_yz_xxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxxz_s_0_0_0[i] += tg_yz_xxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxyy_s_0_0_0[i] += tg_yz_xxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxyz_s_0_0_0[i] += tg_yz_xxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxzz_s_0_0_0[i] += tg_yz_xxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxyyy_s_0_0_0[i] += tg_yz_xxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxyyz_s_0_0_0[i] += tg_yz_xxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxyzz_s_0_0_0[i] += tg_yz_xxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxzzz_s_0_0_0[i] += tg_yz_xxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyyyy_s_0_0_0[i] += tg_yz_xyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyyyz_s_0_0_0[i] += tg_yz_xyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyyzz_s_0_0_0[i] += tg_yz_xyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyzzz_s_0_0_0[i] += tg_yz_xyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xzzzz_s_0_0_0[i] += tg_yz_xzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyyyy_s_0_0_0[i] += tg_yz_yyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyyyz_s_0_0_0[i] += tg_yz_yyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyyzz_s_0_0_0[i] += tg_yz_yyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyzzz_s_0_0_0[i] += tg_yz_yyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yzzzz_s_0_0_0[i] += tg_yz_yzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_zzzzz_s_0_0_0[i] += tg_yz_zzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxxx_s_0_0_0[i] += tg_zz_xxxxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxxy_s_0_0_0[i] += tg_zz_xxxxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxxz_s_0_0_0[i] += tg_zz_xxxxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxyy_s_0_0_0[i] += tg_zz_xxxyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxyz_s_0_0_0[i] += tg_zz_xxxyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxzz_s_0_0_0[i] += tg_zz_xxxzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxyyy_s_0_0_0[i] += tg_zz_xxyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxyyz_s_0_0_0[i] += tg_zz_xxyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxyzz_s_0_0_0[i] += tg_zz_xxyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxzzz_s_0_0_0[i] += tg_zz_xxzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyyyy_s_0_0_0[i] += tg_zz_xyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyyyz_s_0_0_0[i] += tg_zz_xyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyyzz_s_0_0_0[i] += tg_zz_xyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyzzz_s_0_0_0[i] += tg_zz_xyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xzzzz_s_0_0_0[i] += tg_zz_xzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyyyy_s_0_0_0[i] += tg_zz_yyyyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyyyz_s_0_0_0[i] += tg_zz_yyyyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyyzz_s_0_0_0[i] += tg_zz_yyyzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyzzz_s_0_0_0[i] += tg_zz_yyzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yzzzz_s_0_0_0[i] += tg_zz_yzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_zzzzz_s_0_0_0[i] += tg_zz_zzzzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyy_xxxxx_s_0_0_0[i] += tg_y_xxxxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxxy_s_0_0_0[i] += tg_y_xxxxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxxz_s_0_0_0[i] += tg_y_xxxxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxyy_s_0_0_0[i] += tg_y_xxxyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxyz_s_0_0_0[i] += tg_y_xxxyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxzz_s_0_0_0[i] += tg_y_xxxzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxyyy_s_0_0_0[i] += tg_y_xxyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxyyz_s_0_0_0[i] += tg_y_xxyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxyzz_s_0_0_0[i] += tg_y_xxyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxzzz_s_0_0_0[i] += tg_y_xxzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyyyy_s_0_0_0[i] += tg_y_xyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyyyz_s_0_0_0[i] += tg_y_xyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyyzz_s_0_0_0[i] += tg_y_xyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyzzz_s_0_0_0[i] += tg_y_xyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xzzzz_s_0_0_0[i] += tg_y_xzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyyyy_s_0_0_0[i] += tg_y_yyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyyyz_s_0_0_0[i] += tg_y_yyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyyzz_s_0_0_0[i] += tg_y_yyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyzzz_s_0_0_0[i] += tg_y_yyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yzzzz_s_0_0_0[i] += tg_y_yzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_zzzzz_s_0_0_0[i] += tg_y_zzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_zzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyz_xxxxx_s_0_0_0[i] += tg_yy_xxxxx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxxy_s_0_0_0[i] += tg_yy_xxxxy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxxz_s_0_0_0[i] += tg_yy_xxxxz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxyy_s_0_0_0[i] += tg_yy_xxxyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxyz_s_0_0_0[i] += tg_yy_xxxyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxzz_s_0_0_0[i] += tg_yy_xxxzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxyyy_s_0_0_0[i] += tg_yy_xxyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxyyz_s_0_0_0[i] += tg_yy_xxyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxyzz_s_0_0_0[i] += tg_yy_xxyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxzzz_s_0_0_0[i] += tg_yy_xxzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyyyy_s_0_0_0[i] += tg_yy_xyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyyyz_s_0_0_0[i] += tg_yy_xyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyyzz_s_0_0_0[i] += tg_yy_xyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyzzz_s_0_0_0[i] += tg_yy_xyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xzzzz_s_0_0_0[i] += tg_yy_xzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyyyy_s_0_0_0[i] += tg_yy_yyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyyyz_s_0_0_0[i] += tg_yy_yyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyyzz_s_0_0_0[i] += tg_yy_yyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyzzz_s_0_0_0[i] += tg_yy_yyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yzzzz_s_0_0_0[i] += tg_yy_yzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_zzzzz_s_0_0_0[i] += tg_yy_zzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yzz_xxxxx_s_0_0_0[i] += tg_zz_xxxxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxxy_s_0_0_0[i] += tg_zz_xxxxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxxz_s_0_0_0[i] += tg_zz_xxxxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxyy_s_0_0_0[i] += tg_zz_xxxyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxyz_s_0_0_0[i] += tg_zz_xxxyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxzz_s_0_0_0[i] += tg_zz_xxxzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxyyy_s_0_0_0[i] += tg_zz_xxyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxyyz_s_0_0_0[i] += tg_zz_xxyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxyzz_s_0_0_0[i] += tg_zz_xxyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxzzz_s_0_0_0[i] += tg_zz_xxzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyyyy_s_0_0_0[i] += tg_zz_xyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyyyz_s_0_0_0[i] += tg_zz_xyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyyzz_s_0_0_0[i] += tg_zz_xyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyzzz_s_0_0_0[i] += tg_zz_xyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xzzzz_s_0_0_0[i] += tg_zz_xzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyyyy_s_0_0_0[i] += tg_zz_yyyyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyyyz_s_0_0_0[i] += tg_zz_yyyyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyyzz_s_0_0_0[i] += tg_zz_yyyzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyzzz_s_0_0_0[i] += tg_zz_yyzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yzzzz_s_0_0_0[i] += tg_zz_yzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_zzzzz_s_0_0_0[i] += tg_zz_zzzzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzz_xxxxx_s_0_0_0[i] += tg_z_xxxxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxxx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxxy_s_0_0_0[i] += tg_z_xxxxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxxy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxxz_s_0_0_0[i] += tg_z_xxxxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxxz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxyy_s_0_0_0[i] += tg_z_xxxyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxyz_s_0_0_0[i] += tg_z_xxxyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxzz_s_0_0_0[i] += tg_z_xxxzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxyyy_s_0_0_0[i] += tg_z_xxyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxyyz_s_0_0_0[i] += tg_z_xxyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxyzz_s_0_0_0[i] += tg_z_xxyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxzzz_s_0_0_0[i] += tg_z_xxzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyyyy_s_0_0_0[i] += tg_z_xyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyyyz_s_0_0_0[i] += tg_z_xyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyyzz_s_0_0_0[i] += tg_z_xyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyzzz_s_0_0_0[i] += tg_z_xyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xzzzz_s_0_0_0[i] += tg_z_xzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyyyy_s_0_0_0[i] += tg_z_yyyyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyyyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyyyz_s_0_0_0[i] += tg_z_yyyyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyyyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyyzz_s_0_0_0[i] += tg_z_yyyzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyyzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyzzz_s_0_0_0[i] += tg_z_yyzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yzzzz_s_0_0_0[i] += tg_z_yzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_zzzzz_s_0_0_0[i] += tg_z_zzzzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_zzzzz_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace


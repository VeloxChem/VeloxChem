#include "ElectricDipoleMomentumPrimRecIP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_ip(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_ip,
                                      const size_t              idx_dip_gp,
                                      const size_t              idx_dip_hs,
                                      const size_t              idx_ovl_hp,
                                      const size_t              idx_dip_hp,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GP

    auto tr_x_xxxx_x = pbuffer.data(idx_dip_gp);

    auto tr_x_xxxx_y = pbuffer.data(idx_dip_gp + 1);

    auto tr_x_xxxx_z = pbuffer.data(idx_dip_gp + 2);

    auto tr_x_xxxy_x = pbuffer.data(idx_dip_gp + 3);

    auto tr_x_xxxy_z = pbuffer.data(idx_dip_gp + 5);

    auto tr_x_xxxz_x = pbuffer.data(idx_dip_gp + 6);

    auto tr_x_xxxz_y = pbuffer.data(idx_dip_gp + 7);

    auto tr_x_xxxz_z = pbuffer.data(idx_dip_gp + 8);

    auto tr_x_xxyy_x = pbuffer.data(idx_dip_gp + 9);

    auto tr_x_xxyy_y = pbuffer.data(idx_dip_gp + 10);

    auto tr_x_xxyy_z = pbuffer.data(idx_dip_gp + 11);

    auto tr_x_xxyz_z = pbuffer.data(idx_dip_gp + 14);

    auto tr_x_xxzz_x = pbuffer.data(idx_dip_gp + 15);

    auto tr_x_xxzz_y = pbuffer.data(idx_dip_gp + 16);

    auto tr_x_xxzz_z = pbuffer.data(idx_dip_gp + 17);

    auto tr_x_xyyy_x = pbuffer.data(idx_dip_gp + 18);

    auto tr_x_xyyy_y = pbuffer.data(idx_dip_gp + 19);

    auto tr_x_xyzz_x = pbuffer.data(idx_dip_gp + 24);

    auto tr_x_xzzz_x = pbuffer.data(idx_dip_gp + 27);

    auto tr_x_xzzz_z = pbuffer.data(idx_dip_gp + 29);

    auto tr_x_yyyy_x = pbuffer.data(idx_dip_gp + 30);

    auto tr_x_yyyy_y = pbuffer.data(idx_dip_gp + 31);

    auto tr_x_yyyy_z = pbuffer.data(idx_dip_gp + 32);

    auto tr_x_yyyz_y = pbuffer.data(idx_dip_gp + 34);

    auto tr_x_yyyz_z = pbuffer.data(idx_dip_gp + 35);

    auto tr_x_yyzz_x = pbuffer.data(idx_dip_gp + 36);

    auto tr_x_yyzz_y = pbuffer.data(idx_dip_gp + 37);

    auto tr_x_yyzz_z = pbuffer.data(idx_dip_gp + 38);

    auto tr_x_yzzz_x = pbuffer.data(idx_dip_gp + 39);

    auto tr_x_yzzz_z = pbuffer.data(idx_dip_gp + 41);

    auto tr_x_zzzz_x = pbuffer.data(idx_dip_gp + 42);

    auto tr_x_zzzz_y = pbuffer.data(idx_dip_gp + 43);

    auto tr_x_zzzz_z = pbuffer.data(idx_dip_gp + 44);

    auto tr_y_xxxx_x = pbuffer.data(idx_dip_gp + 45);

    auto tr_y_xxxx_y = pbuffer.data(idx_dip_gp + 46);

    auto tr_y_xxxx_z = pbuffer.data(idx_dip_gp + 47);

    auto tr_y_xxxy_y = pbuffer.data(idx_dip_gp + 49);

    auto tr_y_xxxy_z = pbuffer.data(idx_dip_gp + 50);

    auto tr_y_xxxz_x = pbuffer.data(idx_dip_gp + 51);

    auto tr_y_xxxz_z = pbuffer.data(idx_dip_gp + 53);

    auto tr_y_xxyy_x = pbuffer.data(idx_dip_gp + 54);

    auto tr_y_xxyy_y = pbuffer.data(idx_dip_gp + 55);

    auto tr_y_xxyy_z = pbuffer.data(idx_dip_gp + 56);

    auto tr_y_xxyz_z = pbuffer.data(idx_dip_gp + 59);

    auto tr_y_xxzz_x = pbuffer.data(idx_dip_gp + 60);

    auto tr_y_xxzz_y = pbuffer.data(idx_dip_gp + 61);

    auto tr_y_xxzz_z = pbuffer.data(idx_dip_gp + 62);

    auto tr_y_xyyy_x = pbuffer.data(idx_dip_gp + 63);

    auto tr_y_xyyy_y = pbuffer.data(idx_dip_gp + 64);

    auto tr_y_xyyy_z = pbuffer.data(idx_dip_gp + 65);

    auto tr_y_xyyz_z = pbuffer.data(idx_dip_gp + 68);

    auto tr_y_xyzz_y = pbuffer.data(idx_dip_gp + 70);

    auto tr_y_xyzz_z = pbuffer.data(idx_dip_gp + 71);

    auto tr_y_xzzz_y = pbuffer.data(idx_dip_gp + 73);

    auto tr_y_xzzz_z = pbuffer.data(idx_dip_gp + 74);

    auto tr_y_yyyy_x = pbuffer.data(idx_dip_gp + 75);

    auto tr_y_yyyy_y = pbuffer.data(idx_dip_gp + 76);

    auto tr_y_yyyy_z = pbuffer.data(idx_dip_gp + 77);

    auto tr_y_yyyz_x = pbuffer.data(idx_dip_gp + 78);

    auto tr_y_yyyz_y = pbuffer.data(idx_dip_gp + 79);

    auto tr_y_yyyz_z = pbuffer.data(idx_dip_gp + 80);

    auto tr_y_yyzz_x = pbuffer.data(idx_dip_gp + 81);

    auto tr_y_yyzz_y = pbuffer.data(idx_dip_gp + 82);

    auto tr_y_yyzz_z = pbuffer.data(idx_dip_gp + 83);

    auto tr_y_yzzz_y = pbuffer.data(idx_dip_gp + 85);

    auto tr_y_yzzz_z = pbuffer.data(idx_dip_gp + 86);

    auto tr_y_zzzz_x = pbuffer.data(idx_dip_gp + 87);

    auto tr_y_zzzz_y = pbuffer.data(idx_dip_gp + 88);

    auto tr_y_zzzz_z = pbuffer.data(idx_dip_gp + 89);

    auto tr_z_xxxx_x = pbuffer.data(idx_dip_gp + 90);

    auto tr_z_xxxx_y = pbuffer.data(idx_dip_gp + 91);

    auto tr_z_xxxx_z = pbuffer.data(idx_dip_gp + 92);

    auto tr_z_xxxy_x = pbuffer.data(idx_dip_gp + 93);

    auto tr_z_xxxy_y = pbuffer.data(idx_dip_gp + 94);

    auto tr_z_xxxz_x = pbuffer.data(idx_dip_gp + 96);

    auto tr_z_xxxz_y = pbuffer.data(idx_dip_gp + 97);

    auto tr_z_xxxz_z = pbuffer.data(idx_dip_gp + 98);

    auto tr_z_xxyy_x = pbuffer.data(idx_dip_gp + 99);

    auto tr_z_xxyy_y = pbuffer.data(idx_dip_gp + 100);

    auto tr_z_xxyy_z = pbuffer.data(idx_dip_gp + 101);

    auto tr_z_xxyz_x = pbuffer.data(idx_dip_gp + 102);

    auto tr_z_xxyz_y = pbuffer.data(idx_dip_gp + 103);

    auto tr_z_xxzz_x = pbuffer.data(idx_dip_gp + 105);

    auto tr_z_xxzz_y = pbuffer.data(idx_dip_gp + 106);

    auto tr_z_xxzz_z = pbuffer.data(idx_dip_gp + 107);

    auto tr_z_xyyy_y = pbuffer.data(idx_dip_gp + 109);

    auto tr_z_xyyy_z = pbuffer.data(idx_dip_gp + 110);

    auto tr_z_xyyz_y = pbuffer.data(idx_dip_gp + 112);

    auto tr_z_xyyz_z = pbuffer.data(idx_dip_gp + 113);

    auto tr_z_xyzz_y = pbuffer.data(idx_dip_gp + 115);

    auto tr_z_xzzz_x = pbuffer.data(idx_dip_gp + 117);

    auto tr_z_xzzz_y = pbuffer.data(idx_dip_gp + 118);

    auto tr_z_xzzz_z = pbuffer.data(idx_dip_gp + 119);

    auto tr_z_yyyy_x = pbuffer.data(idx_dip_gp + 120);

    auto tr_z_yyyy_y = pbuffer.data(idx_dip_gp + 121);

    auto tr_z_yyyy_z = pbuffer.data(idx_dip_gp + 122);

    auto tr_z_yyyz_x = pbuffer.data(idx_dip_gp + 123);

    auto tr_z_yyyz_y = pbuffer.data(idx_dip_gp + 124);

    auto tr_z_yyyz_z = pbuffer.data(idx_dip_gp + 125);

    auto tr_z_yyzz_x = pbuffer.data(idx_dip_gp + 126);

    auto tr_z_yyzz_y = pbuffer.data(idx_dip_gp + 127);

    auto tr_z_yyzz_z = pbuffer.data(idx_dip_gp + 128);

    auto tr_z_yzzz_x = pbuffer.data(idx_dip_gp + 129);

    auto tr_z_yzzz_y = pbuffer.data(idx_dip_gp + 130);

    auto tr_z_yzzz_z = pbuffer.data(idx_dip_gp + 131);

    auto tr_z_zzzz_x = pbuffer.data(idx_dip_gp + 132);

    auto tr_z_zzzz_y = pbuffer.data(idx_dip_gp + 133);

    auto tr_z_zzzz_z = pbuffer.data(idx_dip_gp + 134);

    // Set up components of auxiliary buffer : HS

    auto tr_x_xxxxx_0 = pbuffer.data(idx_dip_hs);

    auto tr_x_xxxzz_0 = pbuffer.data(idx_dip_hs + 5);

    auto tr_x_xxzzz_0 = pbuffer.data(idx_dip_hs + 9);

    auto tr_x_yyyyy_0 = pbuffer.data(idx_dip_hs + 15);

    auto tr_x_zzzzz_0 = pbuffer.data(idx_dip_hs + 20);

    auto tr_y_xxxxx_0 = pbuffer.data(idx_dip_hs + 21);

    auto tr_y_xxxyy_0 = pbuffer.data(idx_dip_hs + 24);

    auto tr_y_xxyyy_0 = pbuffer.data(idx_dip_hs + 27);

    auto tr_y_xyyyy_0 = pbuffer.data(idx_dip_hs + 31);

    auto tr_y_yyyyy_0 = pbuffer.data(idx_dip_hs + 36);

    auto tr_y_yyyzz_0 = pbuffer.data(idx_dip_hs + 38);

    auto tr_y_yyzzz_0 = pbuffer.data(idx_dip_hs + 39);

    auto tr_y_yzzzz_0 = pbuffer.data(idx_dip_hs + 40);

    auto tr_y_zzzzz_0 = pbuffer.data(idx_dip_hs + 41);

    auto tr_z_xxxxx_0 = pbuffer.data(idx_dip_hs + 42);

    auto tr_z_xxxzz_0 = pbuffer.data(idx_dip_hs + 47);

    auto tr_z_xxzzz_0 = pbuffer.data(idx_dip_hs + 51);

    auto tr_z_xzzzz_0 = pbuffer.data(idx_dip_hs + 56);

    auto tr_z_yyyyy_0 = pbuffer.data(idx_dip_hs + 57);

    auto tr_z_yyyyz_0 = pbuffer.data(idx_dip_hs + 58);

    auto tr_z_yyyzz_0 = pbuffer.data(idx_dip_hs + 59);

    auto tr_z_yyzzz_0 = pbuffer.data(idx_dip_hs + 60);

    auto tr_z_yzzzz_0 = pbuffer.data(idx_dip_hs + 61);

    auto tr_z_zzzzz_0 = pbuffer.data(idx_dip_hs + 62);

    // Set up components of auxiliary buffer : HP

    auto ts_xxxxx_x = pbuffer.data(idx_ovl_hp);

    auto ts_xxxxx_y = pbuffer.data(idx_ovl_hp + 1);

    auto ts_xxxxx_z = pbuffer.data(idx_ovl_hp + 2);

    auto ts_xxxyy_y = pbuffer.data(idx_ovl_hp + 10);

    auto ts_xxxzz_x = pbuffer.data(idx_ovl_hp + 15);

    auto ts_xxxzz_z = pbuffer.data(idx_ovl_hp + 17);

    auto ts_xxyyy_y = pbuffer.data(idx_ovl_hp + 19);

    auto ts_xxzzz_x = pbuffer.data(idx_ovl_hp + 27);

    auto ts_xxzzz_z = pbuffer.data(idx_ovl_hp + 29);

    auto ts_xyyyy_y = pbuffer.data(idx_ovl_hp + 31);

    auto ts_xzzzz_z = pbuffer.data(idx_ovl_hp + 44);

    auto ts_yyyyy_x = pbuffer.data(idx_ovl_hp + 45);

    auto ts_yyyyy_y = pbuffer.data(idx_ovl_hp + 46);

    auto ts_yyyyy_z = pbuffer.data(idx_ovl_hp + 47);

    auto ts_yyyyz_z = pbuffer.data(idx_ovl_hp + 50);

    auto ts_yyyzz_y = pbuffer.data(idx_ovl_hp + 52);

    auto ts_yyyzz_z = pbuffer.data(idx_ovl_hp + 53);

    auto ts_yyzzz_y = pbuffer.data(idx_ovl_hp + 55);

    auto ts_yyzzz_z = pbuffer.data(idx_ovl_hp + 56);

    auto ts_yzzzz_y = pbuffer.data(idx_ovl_hp + 58);

    auto ts_yzzzz_z = pbuffer.data(idx_ovl_hp + 59);

    auto ts_zzzzz_x = pbuffer.data(idx_ovl_hp + 60);

    auto ts_zzzzz_y = pbuffer.data(idx_ovl_hp + 61);

    auto ts_zzzzz_z = pbuffer.data(idx_ovl_hp + 62);

    // Set up components of auxiliary buffer : HP

    auto tr_x_xxxxx_x = pbuffer.data(idx_dip_hp);

    auto tr_x_xxxxx_y = pbuffer.data(idx_dip_hp + 1);

    auto tr_x_xxxxx_z = pbuffer.data(idx_dip_hp + 2);

    auto tr_x_xxxxy_x = pbuffer.data(idx_dip_hp + 3);

    auto tr_x_xxxxy_y = pbuffer.data(idx_dip_hp + 4);

    auto tr_x_xxxxy_z = pbuffer.data(idx_dip_hp + 5);

    auto tr_x_xxxxz_x = pbuffer.data(idx_dip_hp + 6);

    auto tr_x_xxxxz_y = pbuffer.data(idx_dip_hp + 7);

    auto tr_x_xxxxz_z = pbuffer.data(idx_dip_hp + 8);

    auto tr_x_xxxyy_x = pbuffer.data(idx_dip_hp + 9);

    auto tr_x_xxxyy_y = pbuffer.data(idx_dip_hp + 10);

    auto tr_x_xxxyy_z = pbuffer.data(idx_dip_hp + 11);

    auto tr_x_xxxyz_z = pbuffer.data(idx_dip_hp + 14);

    auto tr_x_xxxzz_x = pbuffer.data(idx_dip_hp + 15);

    auto tr_x_xxxzz_y = pbuffer.data(idx_dip_hp + 16);

    auto tr_x_xxxzz_z = pbuffer.data(idx_dip_hp + 17);

    auto tr_x_xxyyy_x = pbuffer.data(idx_dip_hp + 18);

    auto tr_x_xxyyy_y = pbuffer.data(idx_dip_hp + 19);

    auto tr_x_xxyyy_z = pbuffer.data(idx_dip_hp + 20);

    auto tr_x_xxyyz_y = pbuffer.data(idx_dip_hp + 22);

    auto tr_x_xxyyz_z = pbuffer.data(idx_dip_hp + 23);

    auto tr_x_xxyzz_x = pbuffer.data(idx_dip_hp + 24);

    auto tr_x_xxyzz_z = pbuffer.data(idx_dip_hp + 26);

    auto tr_x_xxzzz_x = pbuffer.data(idx_dip_hp + 27);

    auto tr_x_xxzzz_y = pbuffer.data(idx_dip_hp + 28);

    auto tr_x_xxzzz_z = pbuffer.data(idx_dip_hp + 29);

    auto tr_x_xyyyy_x = pbuffer.data(idx_dip_hp + 30);

    auto tr_x_xyyyy_y = pbuffer.data(idx_dip_hp + 31);

    auto tr_x_xyyzz_x = pbuffer.data(idx_dip_hp + 36);

    auto tr_x_xyzzz_x = pbuffer.data(idx_dip_hp + 39);

    auto tr_x_xzzzz_x = pbuffer.data(idx_dip_hp + 42);

    auto tr_x_xzzzz_z = pbuffer.data(idx_dip_hp + 44);

    auto tr_x_yyyyy_x = pbuffer.data(idx_dip_hp + 45);

    auto tr_x_yyyyy_y = pbuffer.data(idx_dip_hp + 46);

    auto tr_x_yyyyy_z = pbuffer.data(idx_dip_hp + 47);

    auto tr_x_yyyyz_y = pbuffer.data(idx_dip_hp + 49);

    auto tr_x_yyyyz_z = pbuffer.data(idx_dip_hp + 50);

    auto tr_x_yyyzz_x = pbuffer.data(idx_dip_hp + 51);

    auto tr_x_yyyzz_y = pbuffer.data(idx_dip_hp + 52);

    auto tr_x_yyyzz_z = pbuffer.data(idx_dip_hp + 53);

    auto tr_x_yyzzz_x = pbuffer.data(idx_dip_hp + 54);

    auto tr_x_yyzzz_y = pbuffer.data(idx_dip_hp + 55);

    auto tr_x_yyzzz_z = pbuffer.data(idx_dip_hp + 56);

    auto tr_x_yzzzz_x = pbuffer.data(idx_dip_hp + 57);

    auto tr_x_yzzzz_y = pbuffer.data(idx_dip_hp + 58);

    auto tr_x_yzzzz_z = pbuffer.data(idx_dip_hp + 59);

    auto tr_x_zzzzz_x = pbuffer.data(idx_dip_hp + 60);

    auto tr_x_zzzzz_y = pbuffer.data(idx_dip_hp + 61);

    auto tr_x_zzzzz_z = pbuffer.data(idx_dip_hp + 62);

    auto tr_y_xxxxx_x = pbuffer.data(idx_dip_hp + 63);

    auto tr_y_xxxxx_y = pbuffer.data(idx_dip_hp + 64);

    auto tr_y_xxxxx_z = pbuffer.data(idx_dip_hp + 65);

    auto tr_y_xxxxy_x = pbuffer.data(idx_dip_hp + 66);

    auto tr_y_xxxxy_y = pbuffer.data(idx_dip_hp + 67);

    auto tr_y_xxxxy_z = pbuffer.data(idx_dip_hp + 68);

    auto tr_y_xxxxz_x = pbuffer.data(idx_dip_hp + 69);

    auto tr_y_xxxxz_z = pbuffer.data(idx_dip_hp + 71);

    auto tr_y_xxxyy_x = pbuffer.data(idx_dip_hp + 72);

    auto tr_y_xxxyy_y = pbuffer.data(idx_dip_hp + 73);

    auto tr_y_xxxyy_z = pbuffer.data(idx_dip_hp + 74);

    auto tr_y_xxxyz_z = pbuffer.data(idx_dip_hp + 77);

    auto tr_y_xxxzz_x = pbuffer.data(idx_dip_hp + 78);

    auto tr_y_xxxzz_y = pbuffer.data(idx_dip_hp + 79);

    auto tr_y_xxxzz_z = pbuffer.data(idx_dip_hp + 80);

    auto tr_y_xxyyy_x = pbuffer.data(idx_dip_hp + 81);

    auto tr_y_xxyyy_y = pbuffer.data(idx_dip_hp + 82);

    auto tr_y_xxyyy_z = pbuffer.data(idx_dip_hp + 83);

    auto tr_y_xxyyz_x = pbuffer.data(idx_dip_hp + 84);

    auto tr_y_xxyyz_z = pbuffer.data(idx_dip_hp + 86);

    auto tr_y_xxyzz_y = pbuffer.data(idx_dip_hp + 88);

    auto tr_y_xxyzz_z = pbuffer.data(idx_dip_hp + 89);

    auto tr_y_xxzzz_x = pbuffer.data(idx_dip_hp + 90);

    auto tr_y_xxzzz_y = pbuffer.data(idx_dip_hp + 91);

    auto tr_y_xxzzz_z = pbuffer.data(idx_dip_hp + 92);

    auto tr_y_xyyyy_x = pbuffer.data(idx_dip_hp + 93);

    auto tr_y_xyyyy_y = pbuffer.data(idx_dip_hp + 94);

    auto tr_y_xyyyy_z = pbuffer.data(idx_dip_hp + 95);

    auto tr_y_xyyyz_z = pbuffer.data(idx_dip_hp + 98);

    auto tr_y_xyyzz_y = pbuffer.data(idx_dip_hp + 100);

    auto tr_y_xyyzz_z = pbuffer.data(idx_dip_hp + 101);

    auto tr_y_xyzzz_y = pbuffer.data(idx_dip_hp + 103);

    auto tr_y_xyzzz_z = pbuffer.data(idx_dip_hp + 104);

    auto tr_y_xzzzz_y = pbuffer.data(idx_dip_hp + 106);

    auto tr_y_xzzzz_z = pbuffer.data(idx_dip_hp + 107);

    auto tr_y_yyyyy_x = pbuffer.data(idx_dip_hp + 108);

    auto tr_y_yyyyy_y = pbuffer.data(idx_dip_hp + 109);

    auto tr_y_yyyyy_z = pbuffer.data(idx_dip_hp + 110);

    auto tr_y_yyyyz_x = pbuffer.data(idx_dip_hp + 111);

    auto tr_y_yyyyz_y = pbuffer.data(idx_dip_hp + 112);

    auto tr_y_yyyyz_z = pbuffer.data(idx_dip_hp + 113);

    auto tr_y_yyyzz_x = pbuffer.data(idx_dip_hp + 114);

    auto tr_y_yyyzz_y = pbuffer.data(idx_dip_hp + 115);

    auto tr_y_yyyzz_z = pbuffer.data(idx_dip_hp + 116);

    auto tr_y_yyzzz_x = pbuffer.data(idx_dip_hp + 117);

    auto tr_y_yyzzz_y = pbuffer.data(idx_dip_hp + 118);

    auto tr_y_yyzzz_z = pbuffer.data(idx_dip_hp + 119);

    auto tr_y_yzzzz_x = pbuffer.data(idx_dip_hp + 120);

    auto tr_y_yzzzz_y = pbuffer.data(idx_dip_hp + 121);

    auto tr_y_yzzzz_z = pbuffer.data(idx_dip_hp + 122);

    auto tr_y_zzzzz_x = pbuffer.data(idx_dip_hp + 123);

    auto tr_y_zzzzz_y = pbuffer.data(idx_dip_hp + 124);

    auto tr_y_zzzzz_z = pbuffer.data(idx_dip_hp + 125);

    auto tr_z_xxxxx_x = pbuffer.data(idx_dip_hp + 126);

    auto tr_z_xxxxx_y = pbuffer.data(idx_dip_hp + 127);

    auto tr_z_xxxxx_z = pbuffer.data(idx_dip_hp + 128);

    auto tr_z_xxxxy_x = pbuffer.data(idx_dip_hp + 129);

    auto tr_z_xxxxy_y = pbuffer.data(idx_dip_hp + 130);

    auto tr_z_xxxxz_x = pbuffer.data(idx_dip_hp + 132);

    auto tr_z_xxxxz_y = pbuffer.data(idx_dip_hp + 133);

    auto tr_z_xxxxz_z = pbuffer.data(idx_dip_hp + 134);

    auto tr_z_xxxyy_x = pbuffer.data(idx_dip_hp + 135);

    auto tr_z_xxxyy_y = pbuffer.data(idx_dip_hp + 136);

    auto tr_z_xxxyy_z = pbuffer.data(idx_dip_hp + 137);

    auto tr_z_xxxyz_x = pbuffer.data(idx_dip_hp + 138);

    auto tr_z_xxxyz_y = pbuffer.data(idx_dip_hp + 139);

    auto tr_z_xxxzz_x = pbuffer.data(idx_dip_hp + 141);

    auto tr_z_xxxzz_y = pbuffer.data(idx_dip_hp + 142);

    auto tr_z_xxxzz_z = pbuffer.data(idx_dip_hp + 143);

    auto tr_z_xxyyy_x = pbuffer.data(idx_dip_hp + 144);

    auto tr_z_xxyyy_y = pbuffer.data(idx_dip_hp + 145);

    auto tr_z_xxyyy_z = pbuffer.data(idx_dip_hp + 146);

    auto tr_z_xxyyz_x = pbuffer.data(idx_dip_hp + 147);

    auto tr_z_xxyyz_y = pbuffer.data(idx_dip_hp + 148);

    auto tr_z_xxyyz_z = pbuffer.data(idx_dip_hp + 149);

    auto tr_z_xxyzz_x = pbuffer.data(idx_dip_hp + 150);

    auto tr_z_xxyzz_y = pbuffer.data(idx_dip_hp + 151);

    auto tr_z_xxzzz_x = pbuffer.data(idx_dip_hp + 153);

    auto tr_z_xxzzz_y = pbuffer.data(idx_dip_hp + 154);

    auto tr_z_xxzzz_z = pbuffer.data(idx_dip_hp + 155);

    auto tr_z_xyyyy_y = pbuffer.data(idx_dip_hp + 157);

    auto tr_z_xyyyy_z = pbuffer.data(idx_dip_hp + 158);

    auto tr_z_xyyyz_y = pbuffer.data(idx_dip_hp + 160);

    auto tr_z_xyyyz_z = pbuffer.data(idx_dip_hp + 161);

    auto tr_z_xyyzz_y = pbuffer.data(idx_dip_hp + 163);

    auto tr_z_xyyzz_z = pbuffer.data(idx_dip_hp + 164);

    auto tr_z_xyzzz_y = pbuffer.data(idx_dip_hp + 166);

    auto tr_z_xzzzz_x = pbuffer.data(idx_dip_hp + 168);

    auto tr_z_xzzzz_y = pbuffer.data(idx_dip_hp + 169);

    auto tr_z_xzzzz_z = pbuffer.data(idx_dip_hp + 170);

    auto tr_z_yyyyy_x = pbuffer.data(idx_dip_hp + 171);

    auto tr_z_yyyyy_y = pbuffer.data(idx_dip_hp + 172);

    auto tr_z_yyyyy_z = pbuffer.data(idx_dip_hp + 173);

    auto tr_z_yyyyz_x = pbuffer.data(idx_dip_hp + 174);

    auto tr_z_yyyyz_y = pbuffer.data(idx_dip_hp + 175);

    auto tr_z_yyyyz_z = pbuffer.data(idx_dip_hp + 176);

    auto tr_z_yyyzz_x = pbuffer.data(idx_dip_hp + 177);

    auto tr_z_yyyzz_y = pbuffer.data(idx_dip_hp + 178);

    auto tr_z_yyyzz_z = pbuffer.data(idx_dip_hp + 179);

    auto tr_z_yyzzz_x = pbuffer.data(idx_dip_hp + 180);

    auto tr_z_yyzzz_y = pbuffer.data(idx_dip_hp + 181);

    auto tr_z_yyzzz_z = pbuffer.data(idx_dip_hp + 182);

    auto tr_z_yzzzz_x = pbuffer.data(idx_dip_hp + 183);

    auto tr_z_yzzzz_y = pbuffer.data(idx_dip_hp + 184);

    auto tr_z_yzzzz_z = pbuffer.data(idx_dip_hp + 185);

    auto tr_z_zzzzz_x = pbuffer.data(idx_dip_hp + 186);

    auto tr_z_zzzzz_y = pbuffer.data(idx_dip_hp + 187);

    auto tr_z_zzzzz_z = pbuffer.data(idx_dip_hp + 188);

    // Set up 0-3 components of targeted buffer : IP

    auto tr_x_xxxxxx_x = pbuffer.data(idx_dip_ip);

    auto tr_x_xxxxxx_y = pbuffer.data(idx_dip_ip + 1);

    auto tr_x_xxxxxx_z = pbuffer.data(idx_dip_ip + 2);

#pragma omp simd aligned(pa_x,              \
                             tr_x_xxxx_x,   \
                             tr_x_xxxx_y,   \
                             tr_x_xxxx_z,   \
                             tr_x_xxxxx_0,  \
                             tr_x_xxxxx_x,  \
                             tr_x_xxxxx_y,  \
                             tr_x_xxxxx_z,  \
                             tr_x_xxxxxx_x, \
                             tr_x_xxxxxx_y, \
                             tr_x_xxxxxx_z, \
                             ts_xxxxx_x,    \
                             ts_xxxxx_y,    \
                             ts_xxxxx_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxx_x[i] = 5.0 * tr_x_xxxx_x[i] * fe_0 + tr_x_xxxxx_0[i] * fe_0 + ts_xxxxx_x[i] * fe_0 + tr_x_xxxxx_x[i] * pa_x[i];

        tr_x_xxxxxx_y[i] = 5.0 * tr_x_xxxx_y[i] * fe_0 + ts_xxxxx_y[i] * fe_0 + tr_x_xxxxx_y[i] * pa_x[i];

        tr_x_xxxxxx_z[i] = 5.0 * tr_x_xxxx_z[i] * fe_0 + ts_xxxxx_z[i] * fe_0 + tr_x_xxxxx_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : IP

    auto tr_x_xxxxxy_x = pbuffer.data(idx_dip_ip + 3);

    auto tr_x_xxxxxy_y = pbuffer.data(idx_dip_ip + 4);

    auto tr_x_xxxxxy_z = pbuffer.data(idx_dip_ip + 5);

#pragma omp simd aligned(pa_y, tr_x_xxxxx_0, tr_x_xxxxx_x, tr_x_xxxxx_y, tr_x_xxxxx_z, tr_x_xxxxxy_x, tr_x_xxxxxy_y, tr_x_xxxxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxy_x[i] = tr_x_xxxxx_x[i] * pa_y[i];

        tr_x_xxxxxy_y[i] = tr_x_xxxxx_0[i] * fe_0 + tr_x_xxxxx_y[i] * pa_y[i];

        tr_x_xxxxxy_z[i] = tr_x_xxxxx_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : IP

    auto tr_x_xxxxxz_x = pbuffer.data(idx_dip_ip + 6);

    auto tr_x_xxxxxz_y = pbuffer.data(idx_dip_ip + 7);

    auto tr_x_xxxxxz_z = pbuffer.data(idx_dip_ip + 8);

#pragma omp simd aligned(pa_z, tr_x_xxxxx_0, tr_x_xxxxx_x, tr_x_xxxxx_y, tr_x_xxxxx_z, tr_x_xxxxxz_x, tr_x_xxxxxz_y, tr_x_xxxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxxz_x[i] = tr_x_xxxxx_x[i] * pa_z[i];

        tr_x_xxxxxz_y[i] = tr_x_xxxxx_y[i] * pa_z[i];

        tr_x_xxxxxz_z[i] = tr_x_xxxxx_0[i] * fe_0 + tr_x_xxxxx_z[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : IP

    auto tr_x_xxxxyy_x = pbuffer.data(idx_dip_ip + 9);

    auto tr_x_xxxxyy_y = pbuffer.data(idx_dip_ip + 10);

    auto tr_x_xxxxyy_z = pbuffer.data(idx_dip_ip + 11);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xxxx_x,   \
                             tr_x_xxxx_z,   \
                             tr_x_xxxxy_x,  \
                             tr_x_xxxxy_z,  \
                             tr_x_xxxxyy_x, \
                             tr_x_xxxxyy_y, \
                             tr_x_xxxxyy_z, \
                             tr_x_xxxyy_y,  \
                             tr_x_xxyy_y,   \
                             ts_xxxyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxyy_x[i] = tr_x_xxxx_x[i] * fe_0 + tr_x_xxxxy_x[i] * pa_y[i];

        tr_x_xxxxyy_y[i] = 3.0 * tr_x_xxyy_y[i] * fe_0 + ts_xxxyy_y[i] * fe_0 + tr_x_xxxyy_y[i] * pa_x[i];

        tr_x_xxxxyy_z[i] = tr_x_xxxx_z[i] * fe_0 + tr_x_xxxxy_z[i] * pa_y[i];
    }

    // Set up 12-15 components of targeted buffer : IP

    auto tr_x_xxxxyz_x = pbuffer.data(idx_dip_ip + 12);

    auto tr_x_xxxxyz_y = pbuffer.data(idx_dip_ip + 13);

    auto tr_x_xxxxyz_z = pbuffer.data(idx_dip_ip + 14);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xxxxy_y, tr_x_xxxxyz_x, tr_x_xxxxyz_y, tr_x_xxxxyz_z, tr_x_xxxxz_x, tr_x_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_x_xxxxyz_x[i] = tr_x_xxxxz_x[i] * pa_y[i];

        tr_x_xxxxyz_y[i] = tr_x_xxxxy_y[i] * pa_z[i];

        tr_x_xxxxyz_z[i] = tr_x_xxxxz_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : IP

    auto tr_x_xxxxzz_x = pbuffer.data(idx_dip_ip + 15);

    auto tr_x_xxxxzz_y = pbuffer.data(idx_dip_ip + 16);

    auto tr_x_xxxxzz_z = pbuffer.data(idx_dip_ip + 17);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_xxxx_x,   \
                             tr_x_xxxx_y,   \
                             tr_x_xxxxz_x,  \
                             tr_x_xxxxz_y,  \
                             tr_x_xxxxzz_x, \
                             tr_x_xxxxzz_y, \
                             tr_x_xxxxzz_z, \
                             tr_x_xxxzz_z,  \
                             tr_x_xxzz_z,   \
                             ts_xxxzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxzz_x[i] = tr_x_xxxx_x[i] * fe_0 + tr_x_xxxxz_x[i] * pa_z[i];

        tr_x_xxxxzz_y[i] = tr_x_xxxx_y[i] * fe_0 + tr_x_xxxxz_y[i] * pa_z[i];

        tr_x_xxxxzz_z[i] = 3.0 * tr_x_xxzz_z[i] * fe_0 + ts_xxxzz_z[i] * fe_0 + tr_x_xxxzz_z[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : IP

    auto tr_x_xxxyyy_x = pbuffer.data(idx_dip_ip + 18);

    auto tr_x_xxxyyy_y = pbuffer.data(idx_dip_ip + 19);

    auto tr_x_xxxyyy_z = pbuffer.data(idx_dip_ip + 20);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xxxy_x,   \
                             tr_x_xxxy_z,   \
                             tr_x_xxxyy_x,  \
                             tr_x_xxxyy_z,  \
                             tr_x_xxxyyy_x, \
                             tr_x_xxxyyy_y, \
                             tr_x_xxxyyy_z, \
                             tr_x_xxyyy_y,  \
                             tr_x_xyyy_y,   \
                             ts_xxyyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyy_x[i] = 2.0 * tr_x_xxxy_x[i] * fe_0 + tr_x_xxxyy_x[i] * pa_y[i];

        tr_x_xxxyyy_y[i] = 2.0 * tr_x_xyyy_y[i] * fe_0 + ts_xxyyy_y[i] * fe_0 + tr_x_xxyyy_y[i] * pa_x[i];

        tr_x_xxxyyy_z[i] = 2.0 * tr_x_xxxy_z[i] * fe_0 + tr_x_xxxyy_z[i] * pa_y[i];
    }

    // Set up 21-24 components of targeted buffer : IP

    auto tr_x_xxxyyz_x = pbuffer.data(idx_dip_ip + 21);

    auto tr_x_xxxyyz_y = pbuffer.data(idx_dip_ip + 22);

    auto tr_x_xxxyyz_z = pbuffer.data(idx_dip_ip + 23);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xxxyy_x, tr_x_xxxyy_y, tr_x_xxxyyz_x, tr_x_xxxyyz_y, tr_x_xxxyyz_z, tr_x_xxxyz_z, tr_x_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyyz_x[i] = tr_x_xxxyy_x[i] * pa_z[i];

        tr_x_xxxyyz_y[i] = tr_x_xxxyy_y[i] * pa_z[i];

        tr_x_xxxyyz_z[i] = tr_x_xxxz_z[i] * fe_0 + tr_x_xxxyz_z[i] * pa_y[i];
    }

    // Set up 24-27 components of targeted buffer : IP

    auto tr_x_xxxyzz_x = pbuffer.data(idx_dip_ip + 24);

    auto tr_x_xxxyzz_y = pbuffer.data(idx_dip_ip + 25);

    auto tr_x_xxxyzz_z = pbuffer.data(idx_dip_ip + 26);

#pragma omp simd aligned(pa_y, tr_x_xxxyzz_x, tr_x_xxxyzz_y, tr_x_xxxyzz_z, tr_x_xxxzz_0, tr_x_xxxzz_x, tr_x_xxxzz_y, tr_x_xxxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyzz_x[i] = tr_x_xxxzz_x[i] * pa_y[i];

        tr_x_xxxyzz_y[i] = tr_x_xxxzz_0[i] * fe_0 + tr_x_xxxzz_y[i] * pa_y[i];

        tr_x_xxxyzz_z[i] = tr_x_xxxzz_z[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : IP

    auto tr_x_xxxzzz_x = pbuffer.data(idx_dip_ip + 27);

    auto tr_x_xxxzzz_y = pbuffer.data(idx_dip_ip + 28);

    auto tr_x_xxxzzz_z = pbuffer.data(idx_dip_ip + 29);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_xxxz_x,   \
                             tr_x_xxxz_y,   \
                             tr_x_xxxzz_x,  \
                             tr_x_xxxzz_y,  \
                             tr_x_xxxzzz_x, \
                             tr_x_xxxzzz_y, \
                             tr_x_xxxzzz_z, \
                             tr_x_xxzzz_z,  \
                             tr_x_xzzz_z,   \
                             ts_xxzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzzz_x[i] = 2.0 * tr_x_xxxz_x[i] * fe_0 + tr_x_xxxzz_x[i] * pa_z[i];

        tr_x_xxxzzz_y[i] = 2.0 * tr_x_xxxz_y[i] * fe_0 + tr_x_xxxzz_y[i] * pa_z[i];

        tr_x_xxxzzz_z[i] = 2.0 * tr_x_xzzz_z[i] * fe_0 + ts_xxzzz_z[i] * fe_0 + tr_x_xxzzz_z[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : IP

    auto tr_x_xxyyyy_x = pbuffer.data(idx_dip_ip + 30);

    auto tr_x_xxyyyy_y = pbuffer.data(idx_dip_ip + 31);

    auto tr_x_xxyyyy_z = pbuffer.data(idx_dip_ip + 32);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xxyy_x,   \
                             tr_x_xxyy_z,   \
                             tr_x_xxyyy_x,  \
                             tr_x_xxyyy_z,  \
                             tr_x_xxyyyy_x, \
                             tr_x_xxyyyy_y, \
                             tr_x_xxyyyy_z, \
                             tr_x_xyyyy_y,  \
                             tr_x_yyyy_y,   \
                             ts_xyyyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyy_x[i] = 3.0 * tr_x_xxyy_x[i] * fe_0 + tr_x_xxyyy_x[i] * pa_y[i];

        tr_x_xxyyyy_y[i] = tr_x_yyyy_y[i] * fe_0 + ts_xyyyy_y[i] * fe_0 + tr_x_xyyyy_y[i] * pa_x[i];

        tr_x_xxyyyy_z[i] = 3.0 * tr_x_xxyy_z[i] * fe_0 + tr_x_xxyyy_z[i] * pa_y[i];
    }

    // Set up 33-36 components of targeted buffer : IP

    auto tr_x_xxyyyz_x = pbuffer.data(idx_dip_ip + 33);

    auto tr_x_xxyyyz_y = pbuffer.data(idx_dip_ip + 34);

    auto tr_x_xxyyyz_z = pbuffer.data(idx_dip_ip + 35);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xxyyy_x, tr_x_xxyyy_y, tr_x_xxyyyz_x, tr_x_xxyyyz_y, tr_x_xxyyyz_z, tr_x_xxyyz_z, tr_x_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyyz_x[i] = tr_x_xxyyy_x[i] * pa_z[i];

        tr_x_xxyyyz_y[i] = tr_x_xxyyy_y[i] * pa_z[i];

        tr_x_xxyyyz_z[i] = 2.0 * tr_x_xxyz_z[i] * fe_0 + tr_x_xxyyz_z[i] * pa_y[i];
    }

    // Set up 36-39 components of targeted buffer : IP

    auto tr_x_xxyyzz_x = pbuffer.data(idx_dip_ip + 36);

    auto tr_x_xxyyzz_y = pbuffer.data(idx_dip_ip + 37);

    auto tr_x_xxyyzz_z = pbuffer.data(idx_dip_ip + 38);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_xxyy_y,   \
                             tr_x_xxyyz_y,  \
                             tr_x_xxyyzz_x, \
                             tr_x_xxyyzz_y, \
                             tr_x_xxyyzz_z, \
                             tr_x_xxyzz_x,  \
                             tr_x_xxyzz_z,  \
                             tr_x_xxzz_x,   \
                             tr_x_xxzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyzz_x[i] = tr_x_xxzz_x[i] * fe_0 + tr_x_xxyzz_x[i] * pa_y[i];

        tr_x_xxyyzz_y[i] = tr_x_xxyy_y[i] * fe_0 + tr_x_xxyyz_y[i] * pa_z[i];

        tr_x_xxyyzz_z[i] = tr_x_xxzz_z[i] * fe_0 + tr_x_xxyzz_z[i] * pa_y[i];
    }

    // Set up 39-42 components of targeted buffer : IP

    auto tr_x_xxyzzz_x = pbuffer.data(idx_dip_ip + 39);

    auto tr_x_xxyzzz_y = pbuffer.data(idx_dip_ip + 40);

    auto tr_x_xxyzzz_z = pbuffer.data(idx_dip_ip + 41);

#pragma omp simd aligned(pa_y, tr_x_xxyzzz_x, tr_x_xxyzzz_y, tr_x_xxyzzz_z, tr_x_xxzzz_0, tr_x_xxzzz_x, tr_x_xxzzz_y, tr_x_xxzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzzz_x[i] = tr_x_xxzzz_x[i] * pa_y[i];

        tr_x_xxyzzz_y[i] = tr_x_xxzzz_0[i] * fe_0 + tr_x_xxzzz_y[i] * pa_y[i];

        tr_x_xxyzzz_z[i] = tr_x_xxzzz_z[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : IP

    auto tr_x_xxzzzz_x = pbuffer.data(idx_dip_ip + 42);

    auto tr_x_xxzzzz_y = pbuffer.data(idx_dip_ip + 43);

    auto tr_x_xxzzzz_z = pbuffer.data(idx_dip_ip + 44);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_xxzz_x,   \
                             tr_x_xxzz_y,   \
                             tr_x_xxzzz_x,  \
                             tr_x_xxzzz_y,  \
                             tr_x_xxzzzz_x, \
                             tr_x_xxzzzz_y, \
                             tr_x_xxzzzz_z, \
                             tr_x_xzzzz_z,  \
                             tr_x_zzzz_z,   \
                             ts_xzzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzzz_x[i] = 3.0 * tr_x_xxzz_x[i] * fe_0 + tr_x_xxzzz_x[i] * pa_z[i];

        tr_x_xxzzzz_y[i] = 3.0 * tr_x_xxzz_y[i] * fe_0 + tr_x_xxzzz_y[i] * pa_z[i];

        tr_x_xxzzzz_z[i] = tr_x_zzzz_z[i] * fe_0 + ts_xzzzz_z[i] * fe_0 + tr_x_xzzzz_z[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : IP

    auto tr_x_xyyyyy_x = pbuffer.data(idx_dip_ip + 45);

    auto tr_x_xyyyyy_y = pbuffer.data(idx_dip_ip + 46);

    auto tr_x_xyyyyy_z = pbuffer.data(idx_dip_ip + 47);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xyyy_x,   \
                             tr_x_xyyyy_x,  \
                             tr_x_xyyyyy_x, \
                             tr_x_xyyyyy_y, \
                             tr_x_xyyyyy_z, \
                             tr_x_yyyyy_y,  \
                             tr_x_yyyyy_z,  \
                             ts_yyyyy_y,    \
                             ts_yyyyy_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyy_x[i] = 4.0 * tr_x_xyyy_x[i] * fe_0 + tr_x_xyyyy_x[i] * pa_y[i];

        tr_x_xyyyyy_y[i] = ts_yyyyy_y[i] * fe_0 + tr_x_yyyyy_y[i] * pa_x[i];

        tr_x_xyyyyy_z[i] = ts_yyyyy_z[i] * fe_0 + tr_x_yyyyy_z[i] * pa_x[i];
    }

    // Set up 48-51 components of targeted buffer : IP

    auto tr_x_xyyyyz_x = pbuffer.data(idx_dip_ip + 48);

    auto tr_x_xyyyyz_y = pbuffer.data(idx_dip_ip + 49);

    auto tr_x_xyyyyz_z = pbuffer.data(idx_dip_ip + 50);

#pragma omp simd aligned(pa_x, pa_z, tr_x_xyyyy_x, tr_x_xyyyy_y, tr_x_xyyyyz_x, tr_x_xyyyyz_y, tr_x_xyyyyz_z, tr_x_yyyyz_z, ts_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyyz_x[i] = tr_x_xyyyy_x[i] * pa_z[i];

        tr_x_xyyyyz_y[i] = tr_x_xyyyy_y[i] * pa_z[i];

        tr_x_xyyyyz_z[i] = ts_yyyyz_z[i] * fe_0 + tr_x_yyyyz_z[i] * pa_x[i];
    }

    // Set up 51-54 components of targeted buffer : IP

    auto tr_x_xyyyzz_x = pbuffer.data(idx_dip_ip + 51);

    auto tr_x_xyyyzz_y = pbuffer.data(idx_dip_ip + 52);

    auto tr_x_xyyyzz_z = pbuffer.data(idx_dip_ip + 53);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xyyyzz_x, \
                             tr_x_xyyyzz_y, \
                             tr_x_xyyyzz_z, \
                             tr_x_xyyzz_x,  \
                             tr_x_xyzz_x,   \
                             tr_x_yyyzz_y,  \
                             tr_x_yyyzz_z,  \
                             ts_yyyzz_y,    \
                             ts_yyyzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyzz_x[i] = 2.0 * tr_x_xyzz_x[i] * fe_0 + tr_x_xyyzz_x[i] * pa_y[i];

        tr_x_xyyyzz_y[i] = ts_yyyzz_y[i] * fe_0 + tr_x_yyyzz_y[i] * pa_x[i];

        tr_x_xyyyzz_z[i] = ts_yyyzz_z[i] * fe_0 + tr_x_yyyzz_z[i] * pa_x[i];
    }

    // Set up 54-57 components of targeted buffer : IP

    auto tr_x_xyyzzz_x = pbuffer.data(idx_dip_ip + 54);

    auto tr_x_xyyzzz_y = pbuffer.data(idx_dip_ip + 55);

    auto tr_x_xyyzzz_z = pbuffer.data(idx_dip_ip + 56);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xyyzzz_x, \
                             tr_x_xyyzzz_y, \
                             tr_x_xyyzzz_z, \
                             tr_x_xyzzz_x,  \
                             tr_x_xzzz_x,   \
                             tr_x_yyzzz_y,  \
                             tr_x_yyzzz_z,  \
                             ts_yyzzz_y,    \
                             ts_yyzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzzz_x[i] = tr_x_xzzz_x[i] * fe_0 + tr_x_xyzzz_x[i] * pa_y[i];

        tr_x_xyyzzz_y[i] = ts_yyzzz_y[i] * fe_0 + tr_x_yyzzz_y[i] * pa_x[i];

        tr_x_xyyzzz_z[i] = ts_yyzzz_z[i] * fe_0 + tr_x_yyzzz_z[i] * pa_x[i];
    }

    // Set up 57-60 components of targeted buffer : IP

    auto tr_x_xyzzzz_x = pbuffer.data(idx_dip_ip + 57);

    auto tr_x_xyzzzz_y = pbuffer.data(idx_dip_ip + 58);

    auto tr_x_xyzzzz_z = pbuffer.data(idx_dip_ip + 59);

#pragma omp simd aligned(pa_x, pa_y, tr_x_xyzzzz_x, tr_x_xyzzzz_y, tr_x_xyzzzz_z, tr_x_xzzzz_x, tr_x_xzzzz_z, tr_x_yzzzz_y, ts_yzzzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzzz_x[i] = tr_x_xzzzz_x[i] * pa_y[i];

        tr_x_xyzzzz_y[i] = ts_yzzzz_y[i] * fe_0 + tr_x_yzzzz_y[i] * pa_x[i];

        tr_x_xyzzzz_z[i] = tr_x_xzzzz_z[i] * pa_y[i];
    }

    // Set up 60-63 components of targeted buffer : IP

    auto tr_x_xzzzzz_x = pbuffer.data(idx_dip_ip + 60);

    auto tr_x_xzzzzz_y = pbuffer.data(idx_dip_ip + 61);

    auto tr_x_xzzzzz_z = pbuffer.data(idx_dip_ip + 62);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_xzzz_x,   \
                             tr_x_xzzzz_x,  \
                             tr_x_xzzzzz_x, \
                             tr_x_xzzzzz_y, \
                             tr_x_xzzzzz_z, \
                             tr_x_zzzzz_y,  \
                             tr_x_zzzzz_z,  \
                             ts_zzzzz_y,    \
                             ts_zzzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzzz_x[i] = 4.0 * tr_x_xzzz_x[i] * fe_0 + tr_x_xzzzz_x[i] * pa_z[i];

        tr_x_xzzzzz_y[i] = ts_zzzzz_y[i] * fe_0 + tr_x_zzzzz_y[i] * pa_x[i];

        tr_x_xzzzzz_z[i] = ts_zzzzz_z[i] * fe_0 + tr_x_zzzzz_z[i] * pa_x[i];
    }

    // Set up 63-66 components of targeted buffer : IP

    auto tr_x_yyyyyy_x = pbuffer.data(idx_dip_ip + 63);

    auto tr_x_yyyyyy_y = pbuffer.data(idx_dip_ip + 64);

    auto tr_x_yyyyyy_z = pbuffer.data(idx_dip_ip + 65);

#pragma omp simd aligned(pa_y,              \
                             tr_x_yyyy_x,   \
                             tr_x_yyyy_y,   \
                             tr_x_yyyy_z,   \
                             tr_x_yyyyy_0,  \
                             tr_x_yyyyy_x,  \
                             tr_x_yyyyy_y,  \
                             tr_x_yyyyy_z,  \
                             tr_x_yyyyyy_x, \
                             tr_x_yyyyyy_y, \
                             tr_x_yyyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyy_x[i] = 5.0 * tr_x_yyyy_x[i] * fe_0 + tr_x_yyyyy_x[i] * pa_y[i];

        tr_x_yyyyyy_y[i] = 5.0 * tr_x_yyyy_y[i] * fe_0 + tr_x_yyyyy_0[i] * fe_0 + tr_x_yyyyy_y[i] * pa_y[i];

        tr_x_yyyyyy_z[i] = 5.0 * tr_x_yyyy_z[i] * fe_0 + tr_x_yyyyy_z[i] * pa_y[i];
    }

    // Set up 66-69 components of targeted buffer : IP

    auto tr_x_yyyyyz_x = pbuffer.data(idx_dip_ip + 66);

    auto tr_x_yyyyyz_y = pbuffer.data(idx_dip_ip + 67);

    auto tr_x_yyyyyz_z = pbuffer.data(idx_dip_ip + 68);

#pragma omp simd aligned(pa_y, pa_z, tr_x_yyyyy_x, tr_x_yyyyy_y, tr_x_yyyyyz_x, tr_x_yyyyyz_y, tr_x_yyyyyz_z, tr_x_yyyyz_z, tr_x_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyyz_x[i] = tr_x_yyyyy_x[i] * pa_z[i];

        tr_x_yyyyyz_y[i] = tr_x_yyyyy_y[i] * pa_z[i];

        tr_x_yyyyyz_z[i] = 4.0 * tr_x_yyyz_z[i] * fe_0 + tr_x_yyyyz_z[i] * pa_y[i];
    }

    // Set up 69-72 components of targeted buffer : IP

    auto tr_x_yyyyzz_x = pbuffer.data(idx_dip_ip + 69);

    auto tr_x_yyyyzz_y = pbuffer.data(idx_dip_ip + 70);

    auto tr_x_yyyyzz_z = pbuffer.data(idx_dip_ip + 71);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_yyyy_y,   \
                             tr_x_yyyyz_y,  \
                             tr_x_yyyyzz_x, \
                             tr_x_yyyyzz_y, \
                             tr_x_yyyyzz_z, \
                             tr_x_yyyzz_x,  \
                             tr_x_yyyzz_z,  \
                             tr_x_yyzz_x,   \
                             tr_x_yyzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyzz_x[i] = 3.0 * tr_x_yyzz_x[i] * fe_0 + tr_x_yyyzz_x[i] * pa_y[i];

        tr_x_yyyyzz_y[i] = tr_x_yyyy_y[i] * fe_0 + tr_x_yyyyz_y[i] * pa_z[i];

        tr_x_yyyyzz_z[i] = 3.0 * tr_x_yyzz_z[i] * fe_0 + tr_x_yyyzz_z[i] * pa_y[i];
    }

    // Set up 72-75 components of targeted buffer : IP

    auto tr_x_yyyzzz_x = pbuffer.data(idx_dip_ip + 72);

    auto tr_x_yyyzzz_y = pbuffer.data(idx_dip_ip + 73);

    auto tr_x_yyyzzz_z = pbuffer.data(idx_dip_ip + 74);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_yyyz_y,   \
                             tr_x_yyyzz_y,  \
                             tr_x_yyyzzz_x, \
                             tr_x_yyyzzz_y, \
                             tr_x_yyyzzz_z, \
                             tr_x_yyzzz_x,  \
                             tr_x_yyzzz_z,  \
                             tr_x_yzzz_x,   \
                             tr_x_yzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzzz_x[i] = 2.0 * tr_x_yzzz_x[i] * fe_0 + tr_x_yyzzz_x[i] * pa_y[i];

        tr_x_yyyzzz_y[i] = 2.0 * tr_x_yyyz_y[i] * fe_0 + tr_x_yyyzz_y[i] * pa_z[i];

        tr_x_yyyzzz_z[i] = 2.0 * tr_x_yzzz_z[i] * fe_0 + tr_x_yyzzz_z[i] * pa_y[i];
    }

    // Set up 75-78 components of targeted buffer : IP

    auto tr_x_yyzzzz_x = pbuffer.data(idx_dip_ip + 75);

    auto tr_x_yyzzzz_y = pbuffer.data(idx_dip_ip + 76);

    auto tr_x_yyzzzz_z = pbuffer.data(idx_dip_ip + 77);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_yyzz_y,   \
                             tr_x_yyzzz_y,  \
                             tr_x_yyzzzz_x, \
                             tr_x_yyzzzz_y, \
                             tr_x_yyzzzz_z, \
                             tr_x_yzzzz_x,  \
                             tr_x_yzzzz_z,  \
                             tr_x_zzzz_x,   \
                             tr_x_zzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzzz_x[i] = tr_x_zzzz_x[i] * fe_0 + tr_x_yzzzz_x[i] * pa_y[i];

        tr_x_yyzzzz_y[i] = 3.0 * tr_x_yyzz_y[i] * fe_0 + tr_x_yyzzz_y[i] * pa_z[i];

        tr_x_yyzzzz_z[i] = tr_x_zzzz_z[i] * fe_0 + tr_x_yzzzz_z[i] * pa_y[i];
    }

    // Set up 78-81 components of targeted buffer : IP

    auto tr_x_yzzzzz_x = pbuffer.data(idx_dip_ip + 78);

    auto tr_x_yzzzzz_y = pbuffer.data(idx_dip_ip + 79);

    auto tr_x_yzzzzz_z = pbuffer.data(idx_dip_ip + 80);

#pragma omp simd aligned(pa_y, tr_x_yzzzzz_x, tr_x_yzzzzz_y, tr_x_yzzzzz_z, tr_x_zzzzz_0, tr_x_zzzzz_x, tr_x_zzzzz_y, tr_x_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzzz_x[i] = tr_x_zzzzz_x[i] * pa_y[i];

        tr_x_yzzzzz_y[i] = tr_x_zzzzz_0[i] * fe_0 + tr_x_zzzzz_y[i] * pa_y[i];

        tr_x_yzzzzz_z[i] = tr_x_zzzzz_z[i] * pa_y[i];
    }

    // Set up 81-84 components of targeted buffer : IP

    auto tr_x_zzzzzz_x = pbuffer.data(idx_dip_ip + 81);

    auto tr_x_zzzzzz_y = pbuffer.data(idx_dip_ip + 82);

    auto tr_x_zzzzzz_z = pbuffer.data(idx_dip_ip + 83);

#pragma omp simd aligned(pa_z,              \
                             tr_x_zzzz_x,   \
                             tr_x_zzzz_y,   \
                             tr_x_zzzz_z,   \
                             tr_x_zzzzz_0,  \
                             tr_x_zzzzz_x,  \
                             tr_x_zzzzz_y,  \
                             tr_x_zzzzz_z,  \
                             tr_x_zzzzzz_x, \
                             tr_x_zzzzzz_y, \
                             tr_x_zzzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzzz_x[i] = 5.0 * tr_x_zzzz_x[i] * fe_0 + tr_x_zzzzz_x[i] * pa_z[i];

        tr_x_zzzzzz_y[i] = 5.0 * tr_x_zzzz_y[i] * fe_0 + tr_x_zzzzz_y[i] * pa_z[i];

        tr_x_zzzzzz_z[i] = 5.0 * tr_x_zzzz_z[i] * fe_0 + tr_x_zzzzz_0[i] * fe_0 + tr_x_zzzzz_z[i] * pa_z[i];
    }

    // Set up 84-87 components of targeted buffer : IP

    auto tr_y_xxxxxx_x = pbuffer.data(idx_dip_ip + 84);

    auto tr_y_xxxxxx_y = pbuffer.data(idx_dip_ip + 85);

    auto tr_y_xxxxxx_z = pbuffer.data(idx_dip_ip + 86);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xxxx_x,   \
                             tr_y_xxxx_y,   \
                             tr_y_xxxx_z,   \
                             tr_y_xxxxx_0,  \
                             tr_y_xxxxx_x,  \
                             tr_y_xxxxx_y,  \
                             tr_y_xxxxx_z,  \
                             tr_y_xxxxxx_x, \
                             tr_y_xxxxxx_y, \
                             tr_y_xxxxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxx_x[i] = 5.0 * tr_y_xxxx_x[i] * fe_0 + tr_y_xxxxx_0[i] * fe_0 + tr_y_xxxxx_x[i] * pa_x[i];

        tr_y_xxxxxx_y[i] = 5.0 * tr_y_xxxx_y[i] * fe_0 + tr_y_xxxxx_y[i] * pa_x[i];

        tr_y_xxxxxx_z[i] = 5.0 * tr_y_xxxx_z[i] * fe_0 + tr_y_xxxxx_z[i] * pa_x[i];
    }

    // Set up 87-90 components of targeted buffer : IP

    auto tr_y_xxxxxy_x = pbuffer.data(idx_dip_ip + 87);

    auto tr_y_xxxxxy_y = pbuffer.data(idx_dip_ip + 88);

    auto tr_y_xxxxxy_z = pbuffer.data(idx_dip_ip + 89);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_y_xxxxx_x,  \
                             tr_y_xxxxxy_x, \
                             tr_y_xxxxxy_y, \
                             tr_y_xxxxxy_z, \
                             tr_y_xxxxy_y,  \
                             tr_y_xxxxy_z,  \
                             tr_y_xxxy_y,   \
                             tr_y_xxxy_z,   \
                             ts_xxxxx_x,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxy_x[i] = ts_xxxxx_x[i] * fe_0 + tr_y_xxxxx_x[i] * pa_y[i];

        tr_y_xxxxxy_y[i] = 4.0 * tr_y_xxxy_y[i] * fe_0 + tr_y_xxxxy_y[i] * pa_x[i];

        tr_y_xxxxxy_z[i] = 4.0 * tr_y_xxxy_z[i] * fe_0 + tr_y_xxxxy_z[i] * pa_x[i];
    }

    // Set up 90-93 components of targeted buffer : IP

    auto tr_y_xxxxxz_x = pbuffer.data(idx_dip_ip + 90);

    auto tr_y_xxxxxz_y = pbuffer.data(idx_dip_ip + 91);

    auto tr_y_xxxxxz_z = pbuffer.data(idx_dip_ip + 92);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxxxx_x, tr_y_xxxxx_y, tr_y_xxxxxz_x, tr_y_xxxxxz_y, tr_y_xxxxxz_z, tr_y_xxxxz_z, tr_y_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxxz_x[i] = tr_y_xxxxx_x[i] * pa_z[i];

        tr_y_xxxxxz_y[i] = tr_y_xxxxx_y[i] * pa_z[i];

        tr_y_xxxxxz_z[i] = 4.0 * tr_y_xxxz_z[i] * fe_0 + tr_y_xxxxz_z[i] * pa_x[i];
    }

    // Set up 93-96 components of targeted buffer : IP

    auto tr_y_xxxxyy_x = pbuffer.data(idx_dip_ip + 93);

    auto tr_y_xxxxyy_y = pbuffer.data(idx_dip_ip + 94);

    auto tr_y_xxxxyy_z = pbuffer.data(idx_dip_ip + 95);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xxxxyy_x, \
                             tr_y_xxxxyy_y, \
                             tr_y_xxxxyy_z, \
                             tr_y_xxxyy_0,  \
                             tr_y_xxxyy_x,  \
                             tr_y_xxxyy_y,  \
                             tr_y_xxxyy_z,  \
                             tr_y_xxyy_x,   \
                             tr_y_xxyy_y,   \
                             tr_y_xxyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyy_x[i] = 3.0 * tr_y_xxyy_x[i] * fe_0 + tr_y_xxxyy_0[i] * fe_0 + tr_y_xxxyy_x[i] * pa_x[i];

        tr_y_xxxxyy_y[i] = 3.0 * tr_y_xxyy_y[i] * fe_0 + tr_y_xxxyy_y[i] * pa_x[i];

        tr_y_xxxxyy_z[i] = 3.0 * tr_y_xxyy_z[i] * fe_0 + tr_y_xxxyy_z[i] * pa_x[i];
    }

    // Set up 96-99 components of targeted buffer : IP

    auto tr_y_xxxxyz_x = pbuffer.data(idx_dip_ip + 96);

    auto tr_y_xxxxyz_y = pbuffer.data(idx_dip_ip + 97);

    auto tr_y_xxxxyz_z = pbuffer.data(idx_dip_ip + 98);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxxxy_x, tr_y_xxxxy_y, tr_y_xxxxyz_x, tr_y_xxxxyz_y, tr_y_xxxxyz_z, tr_y_xxxyz_z, tr_y_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxyz_x[i] = tr_y_xxxxy_x[i] * pa_z[i];

        tr_y_xxxxyz_y[i] = tr_y_xxxxy_y[i] * pa_z[i];

        tr_y_xxxxyz_z[i] = 3.0 * tr_y_xxyz_z[i] * fe_0 + tr_y_xxxyz_z[i] * pa_x[i];
    }

    // Set up 99-102 components of targeted buffer : IP

    auto tr_y_xxxxzz_x = pbuffer.data(idx_dip_ip + 99);

    auto tr_y_xxxxzz_y = pbuffer.data(idx_dip_ip + 100);

    auto tr_y_xxxxzz_z = pbuffer.data(idx_dip_ip + 101);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xxxx_x,   \
                             tr_y_xxxxz_x,  \
                             tr_y_xxxxzz_x, \
                             tr_y_xxxxzz_y, \
                             tr_y_xxxxzz_z, \
                             tr_y_xxxzz_y,  \
                             tr_y_xxxzz_z,  \
                             tr_y_xxzz_y,   \
                             tr_y_xxzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxzz_x[i] = tr_y_xxxx_x[i] * fe_0 + tr_y_xxxxz_x[i] * pa_z[i];

        tr_y_xxxxzz_y[i] = 3.0 * tr_y_xxzz_y[i] * fe_0 + tr_y_xxxzz_y[i] * pa_x[i];

        tr_y_xxxxzz_z[i] = 3.0 * tr_y_xxzz_z[i] * fe_0 + tr_y_xxxzz_z[i] * pa_x[i];
    }

    // Set up 102-105 components of targeted buffer : IP

    auto tr_y_xxxyyy_x = pbuffer.data(idx_dip_ip + 102);

    auto tr_y_xxxyyy_y = pbuffer.data(idx_dip_ip + 103);

    auto tr_y_xxxyyy_z = pbuffer.data(idx_dip_ip + 104);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xxxyyy_x, \
                             tr_y_xxxyyy_y, \
                             tr_y_xxxyyy_z, \
                             tr_y_xxyyy_0,  \
                             tr_y_xxyyy_x,  \
                             tr_y_xxyyy_y,  \
                             tr_y_xxyyy_z,  \
                             tr_y_xyyy_x,   \
                             tr_y_xyyy_y,   \
                             tr_y_xyyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyy_x[i] = 2.0 * tr_y_xyyy_x[i] * fe_0 + tr_y_xxyyy_0[i] * fe_0 + tr_y_xxyyy_x[i] * pa_x[i];

        tr_y_xxxyyy_y[i] = 2.0 * tr_y_xyyy_y[i] * fe_0 + tr_y_xxyyy_y[i] * pa_x[i];

        tr_y_xxxyyy_z[i] = 2.0 * tr_y_xyyy_z[i] * fe_0 + tr_y_xxyyy_z[i] * pa_x[i];
    }

    // Set up 105-108 components of targeted buffer : IP

    auto tr_y_xxxyyz_x = pbuffer.data(idx_dip_ip + 105);

    auto tr_y_xxxyyz_y = pbuffer.data(idx_dip_ip + 106);

    auto tr_y_xxxyyz_z = pbuffer.data(idx_dip_ip + 107);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxxyy_x, tr_y_xxxyy_y, tr_y_xxxyyz_x, tr_y_xxxyyz_y, tr_y_xxxyyz_z, tr_y_xxyyz_z, tr_y_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyyz_x[i] = tr_y_xxxyy_x[i] * pa_z[i];

        tr_y_xxxyyz_y[i] = tr_y_xxxyy_y[i] * pa_z[i];

        tr_y_xxxyyz_z[i] = 2.0 * tr_y_xyyz_z[i] * fe_0 + tr_y_xxyyz_z[i] * pa_x[i];
    }

    // Set up 108-111 components of targeted buffer : IP

    auto tr_y_xxxyzz_x = pbuffer.data(idx_dip_ip + 108);

    auto tr_y_xxxyzz_y = pbuffer.data(idx_dip_ip + 109);

    auto tr_y_xxxyzz_z = pbuffer.data(idx_dip_ip + 110);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_y_xxxyzz_x, \
                             tr_y_xxxyzz_y, \
                             tr_y_xxxyzz_z, \
                             tr_y_xxxzz_x,  \
                             tr_y_xxyzz_y,  \
                             tr_y_xxyzz_z,  \
                             tr_y_xyzz_y,   \
                             tr_y_xyzz_z,   \
                             ts_xxxzz_x,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyzz_x[i] = ts_xxxzz_x[i] * fe_0 + tr_y_xxxzz_x[i] * pa_y[i];

        tr_y_xxxyzz_y[i] = 2.0 * tr_y_xyzz_y[i] * fe_0 + tr_y_xxyzz_y[i] * pa_x[i];

        tr_y_xxxyzz_z[i] = 2.0 * tr_y_xyzz_z[i] * fe_0 + tr_y_xxyzz_z[i] * pa_x[i];
    }

    // Set up 111-114 components of targeted buffer : IP

    auto tr_y_xxxzzz_x = pbuffer.data(idx_dip_ip + 111);

    auto tr_y_xxxzzz_y = pbuffer.data(idx_dip_ip + 112);

    auto tr_y_xxxzzz_z = pbuffer.data(idx_dip_ip + 113);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xxxz_x,   \
                             tr_y_xxxzz_x,  \
                             tr_y_xxxzzz_x, \
                             tr_y_xxxzzz_y, \
                             tr_y_xxxzzz_z, \
                             tr_y_xxzzz_y,  \
                             tr_y_xxzzz_z,  \
                             tr_y_xzzz_y,   \
                             tr_y_xzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzzz_x[i] = 2.0 * tr_y_xxxz_x[i] * fe_0 + tr_y_xxxzz_x[i] * pa_z[i];

        tr_y_xxxzzz_y[i] = 2.0 * tr_y_xzzz_y[i] * fe_0 + tr_y_xxzzz_y[i] * pa_x[i];

        tr_y_xxxzzz_z[i] = 2.0 * tr_y_xzzz_z[i] * fe_0 + tr_y_xxzzz_z[i] * pa_x[i];
    }

    // Set up 114-117 components of targeted buffer : IP

    auto tr_y_xxyyyy_x = pbuffer.data(idx_dip_ip + 114);

    auto tr_y_xxyyyy_y = pbuffer.data(idx_dip_ip + 115);

    auto tr_y_xxyyyy_z = pbuffer.data(idx_dip_ip + 116);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xxyyyy_x, \
                             tr_y_xxyyyy_y, \
                             tr_y_xxyyyy_z, \
                             tr_y_xyyyy_0,  \
                             tr_y_xyyyy_x,  \
                             tr_y_xyyyy_y,  \
                             tr_y_xyyyy_z,  \
                             tr_y_yyyy_x,   \
                             tr_y_yyyy_y,   \
                             tr_y_yyyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyy_x[i] = tr_y_yyyy_x[i] * fe_0 + tr_y_xyyyy_0[i] * fe_0 + tr_y_xyyyy_x[i] * pa_x[i];

        tr_y_xxyyyy_y[i] = tr_y_yyyy_y[i] * fe_0 + tr_y_xyyyy_y[i] * pa_x[i];

        tr_y_xxyyyy_z[i] = tr_y_yyyy_z[i] * fe_0 + tr_y_xyyyy_z[i] * pa_x[i];
    }

    // Set up 117-120 components of targeted buffer : IP

    auto tr_y_xxyyyz_x = pbuffer.data(idx_dip_ip + 117);

    auto tr_y_xxyyyz_y = pbuffer.data(idx_dip_ip + 118);

    auto tr_y_xxyyyz_z = pbuffer.data(idx_dip_ip + 119);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxyyy_x, tr_y_xxyyy_y, tr_y_xxyyyz_x, tr_y_xxyyyz_y, tr_y_xxyyyz_z, tr_y_xyyyz_z, tr_y_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyyz_x[i] = tr_y_xxyyy_x[i] * pa_z[i];

        tr_y_xxyyyz_y[i] = tr_y_xxyyy_y[i] * pa_z[i];

        tr_y_xxyyyz_z[i] = tr_y_yyyz_z[i] * fe_0 + tr_y_xyyyz_z[i] * pa_x[i];
    }

    // Set up 120-123 components of targeted buffer : IP

    auto tr_y_xxyyzz_x = pbuffer.data(idx_dip_ip + 120);

    auto tr_y_xxyyzz_y = pbuffer.data(idx_dip_ip + 121);

    auto tr_y_xxyyzz_z = pbuffer.data(idx_dip_ip + 122);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xxyy_x,   \
                             tr_y_xxyyz_x,  \
                             tr_y_xxyyzz_x, \
                             tr_y_xxyyzz_y, \
                             tr_y_xxyyzz_z, \
                             tr_y_xyyzz_y,  \
                             tr_y_xyyzz_z,  \
                             tr_y_yyzz_y,   \
                             tr_y_yyzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyzz_x[i] = tr_y_xxyy_x[i] * fe_0 + tr_y_xxyyz_x[i] * pa_z[i];

        tr_y_xxyyzz_y[i] = tr_y_yyzz_y[i] * fe_0 + tr_y_xyyzz_y[i] * pa_x[i];

        tr_y_xxyyzz_z[i] = tr_y_yyzz_z[i] * fe_0 + tr_y_xyyzz_z[i] * pa_x[i];
    }

    // Set up 123-126 components of targeted buffer : IP

    auto tr_y_xxyzzz_x = pbuffer.data(idx_dip_ip + 123);

    auto tr_y_xxyzzz_y = pbuffer.data(idx_dip_ip + 124);

    auto tr_y_xxyzzz_z = pbuffer.data(idx_dip_ip + 125);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_y_xxyzzz_x, \
                             tr_y_xxyzzz_y, \
                             tr_y_xxyzzz_z, \
                             tr_y_xxzzz_x,  \
                             tr_y_xyzzz_y,  \
                             tr_y_xyzzz_z,  \
                             tr_y_yzzz_y,   \
                             tr_y_yzzz_z,   \
                             ts_xxzzz_x,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzzz_x[i] = ts_xxzzz_x[i] * fe_0 + tr_y_xxzzz_x[i] * pa_y[i];

        tr_y_xxyzzz_y[i] = tr_y_yzzz_y[i] * fe_0 + tr_y_xyzzz_y[i] * pa_x[i];

        tr_y_xxyzzz_z[i] = tr_y_yzzz_z[i] * fe_0 + tr_y_xyzzz_z[i] * pa_x[i];
    }

    // Set up 126-129 components of targeted buffer : IP

    auto tr_y_xxzzzz_x = pbuffer.data(idx_dip_ip + 126);

    auto tr_y_xxzzzz_y = pbuffer.data(idx_dip_ip + 127);

    auto tr_y_xxzzzz_z = pbuffer.data(idx_dip_ip + 128);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xxzz_x,   \
                             tr_y_xxzzz_x,  \
                             tr_y_xxzzzz_x, \
                             tr_y_xxzzzz_y, \
                             tr_y_xxzzzz_z, \
                             tr_y_xzzzz_y,  \
                             tr_y_xzzzz_z,  \
                             tr_y_zzzz_y,   \
                             tr_y_zzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzzz_x[i] = 3.0 * tr_y_xxzz_x[i] * fe_0 + tr_y_xxzzz_x[i] * pa_z[i];

        tr_y_xxzzzz_y[i] = tr_y_zzzz_y[i] * fe_0 + tr_y_xzzzz_y[i] * pa_x[i];

        tr_y_xxzzzz_z[i] = tr_y_zzzz_z[i] * fe_0 + tr_y_xzzzz_z[i] * pa_x[i];
    }

    // Set up 129-132 components of targeted buffer : IP

    auto tr_y_xyyyyy_x = pbuffer.data(idx_dip_ip + 129);

    auto tr_y_xyyyyy_y = pbuffer.data(idx_dip_ip + 130);

    auto tr_y_xyyyyy_z = pbuffer.data(idx_dip_ip + 131);

#pragma omp simd aligned(pa_x, tr_y_xyyyyy_x, tr_y_xyyyyy_y, tr_y_xyyyyy_z, tr_y_yyyyy_0, tr_y_yyyyy_x, tr_y_yyyyy_y, tr_y_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyyy_x[i] = tr_y_yyyyy_0[i] * fe_0 + tr_y_yyyyy_x[i] * pa_x[i];

        tr_y_xyyyyy_y[i] = tr_y_yyyyy_y[i] * pa_x[i];

        tr_y_xyyyyy_z[i] = tr_y_yyyyy_z[i] * pa_x[i];
    }

    // Set up 132-135 components of targeted buffer : IP

    auto tr_y_xyyyyz_x = pbuffer.data(idx_dip_ip + 132);

    auto tr_y_xyyyyz_y = pbuffer.data(idx_dip_ip + 133);

    auto tr_y_xyyyyz_z = pbuffer.data(idx_dip_ip + 134);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xyyyy_x, tr_y_xyyyyz_x, tr_y_xyyyyz_y, tr_y_xyyyyz_z, tr_y_yyyyz_y, tr_y_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_y_xyyyyz_x[i] = tr_y_xyyyy_x[i] * pa_z[i];

        tr_y_xyyyyz_y[i] = tr_y_yyyyz_y[i] * pa_x[i];

        tr_y_xyyyyz_z[i] = tr_y_yyyyz_z[i] * pa_x[i];
    }

    // Set up 135-138 components of targeted buffer : IP

    auto tr_y_xyyyzz_x = pbuffer.data(idx_dip_ip + 135);

    auto tr_y_xyyyzz_y = pbuffer.data(idx_dip_ip + 136);

    auto tr_y_xyyyzz_z = pbuffer.data(idx_dip_ip + 137);

#pragma omp simd aligned(pa_x, tr_y_xyyyzz_x, tr_y_xyyyzz_y, tr_y_xyyyzz_z, tr_y_yyyzz_0, tr_y_yyyzz_x, tr_y_yyyzz_y, tr_y_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyzz_x[i] = tr_y_yyyzz_0[i] * fe_0 + tr_y_yyyzz_x[i] * pa_x[i];

        tr_y_xyyyzz_y[i] = tr_y_yyyzz_y[i] * pa_x[i];

        tr_y_xyyyzz_z[i] = tr_y_yyyzz_z[i] * pa_x[i];
    }

    // Set up 138-141 components of targeted buffer : IP

    auto tr_y_xyyzzz_x = pbuffer.data(idx_dip_ip + 138);

    auto tr_y_xyyzzz_y = pbuffer.data(idx_dip_ip + 139);

    auto tr_y_xyyzzz_z = pbuffer.data(idx_dip_ip + 140);

#pragma omp simd aligned(pa_x, tr_y_xyyzzz_x, tr_y_xyyzzz_y, tr_y_xyyzzz_z, tr_y_yyzzz_0, tr_y_yyzzz_x, tr_y_yyzzz_y, tr_y_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzzz_x[i] = tr_y_yyzzz_0[i] * fe_0 + tr_y_yyzzz_x[i] * pa_x[i];

        tr_y_xyyzzz_y[i] = tr_y_yyzzz_y[i] * pa_x[i];

        tr_y_xyyzzz_z[i] = tr_y_yyzzz_z[i] * pa_x[i];
    }

    // Set up 141-144 components of targeted buffer : IP

    auto tr_y_xyzzzz_x = pbuffer.data(idx_dip_ip + 141);

    auto tr_y_xyzzzz_y = pbuffer.data(idx_dip_ip + 142);

    auto tr_y_xyzzzz_z = pbuffer.data(idx_dip_ip + 143);

#pragma omp simd aligned(pa_x, tr_y_xyzzzz_x, tr_y_xyzzzz_y, tr_y_xyzzzz_z, tr_y_yzzzz_0, tr_y_yzzzz_x, tr_y_yzzzz_y, tr_y_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzzz_x[i] = tr_y_yzzzz_0[i] * fe_0 + tr_y_yzzzz_x[i] * pa_x[i];

        tr_y_xyzzzz_y[i] = tr_y_yzzzz_y[i] * pa_x[i];

        tr_y_xyzzzz_z[i] = tr_y_yzzzz_z[i] * pa_x[i];
    }

    // Set up 144-147 components of targeted buffer : IP

    auto tr_y_xzzzzz_x = pbuffer.data(idx_dip_ip + 144);

    auto tr_y_xzzzzz_y = pbuffer.data(idx_dip_ip + 145);

    auto tr_y_xzzzzz_z = pbuffer.data(idx_dip_ip + 146);

#pragma omp simd aligned(pa_x, tr_y_xzzzzz_x, tr_y_xzzzzz_y, tr_y_xzzzzz_z, tr_y_zzzzz_0, tr_y_zzzzz_x, tr_y_zzzzz_y, tr_y_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzzz_x[i] = tr_y_zzzzz_0[i] * fe_0 + tr_y_zzzzz_x[i] * pa_x[i];

        tr_y_xzzzzz_y[i] = tr_y_zzzzz_y[i] * pa_x[i];

        tr_y_xzzzzz_z[i] = tr_y_zzzzz_z[i] * pa_x[i];
    }

    // Set up 147-150 components of targeted buffer : IP

    auto tr_y_yyyyyy_x = pbuffer.data(idx_dip_ip + 147);

    auto tr_y_yyyyyy_y = pbuffer.data(idx_dip_ip + 148);

    auto tr_y_yyyyyy_z = pbuffer.data(idx_dip_ip + 149);

#pragma omp simd aligned(pa_y,              \
                             tr_y_yyyy_x,   \
                             tr_y_yyyy_y,   \
                             tr_y_yyyy_z,   \
                             tr_y_yyyyy_0,  \
                             tr_y_yyyyy_x,  \
                             tr_y_yyyyy_y,  \
                             tr_y_yyyyy_z,  \
                             tr_y_yyyyyy_x, \
                             tr_y_yyyyyy_y, \
                             tr_y_yyyyyy_z, \
                             ts_yyyyy_x,    \
                             ts_yyyyy_y,    \
                             ts_yyyyy_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyy_x[i] = 5.0 * tr_y_yyyy_x[i] * fe_0 + ts_yyyyy_x[i] * fe_0 + tr_y_yyyyy_x[i] * pa_y[i];

        tr_y_yyyyyy_y[i] = 5.0 * tr_y_yyyy_y[i] * fe_0 + tr_y_yyyyy_0[i] * fe_0 + ts_yyyyy_y[i] * fe_0 + tr_y_yyyyy_y[i] * pa_y[i];

        tr_y_yyyyyy_z[i] = 5.0 * tr_y_yyyy_z[i] * fe_0 + ts_yyyyy_z[i] * fe_0 + tr_y_yyyyy_z[i] * pa_y[i];
    }

    // Set up 150-153 components of targeted buffer : IP

    auto tr_y_yyyyyz_x = pbuffer.data(idx_dip_ip + 150);

    auto tr_y_yyyyyz_y = pbuffer.data(idx_dip_ip + 151);

    auto tr_y_yyyyyz_z = pbuffer.data(idx_dip_ip + 152);

#pragma omp simd aligned(pa_z, tr_y_yyyyy_0, tr_y_yyyyy_x, tr_y_yyyyy_y, tr_y_yyyyy_z, tr_y_yyyyyz_x, tr_y_yyyyyz_y, tr_y_yyyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyyz_x[i] = tr_y_yyyyy_x[i] * pa_z[i];

        tr_y_yyyyyz_y[i] = tr_y_yyyyy_y[i] * pa_z[i];

        tr_y_yyyyyz_z[i] = tr_y_yyyyy_0[i] * fe_0 + tr_y_yyyyy_z[i] * pa_z[i];
    }

    // Set up 153-156 components of targeted buffer : IP

    auto tr_y_yyyyzz_x = pbuffer.data(idx_dip_ip + 153);

    auto tr_y_yyyyzz_y = pbuffer.data(idx_dip_ip + 154);

    auto tr_y_yyyyzz_z = pbuffer.data(idx_dip_ip + 155);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_yyyy_x,   \
                             tr_y_yyyy_y,   \
                             tr_y_yyyyz_x,  \
                             tr_y_yyyyz_y,  \
                             tr_y_yyyyzz_x, \
                             tr_y_yyyyzz_y, \
                             tr_y_yyyyzz_z, \
                             tr_y_yyyzz_z,  \
                             tr_y_yyzz_z,   \
                             ts_yyyzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyzz_x[i] = tr_y_yyyy_x[i] * fe_0 + tr_y_yyyyz_x[i] * pa_z[i];

        tr_y_yyyyzz_y[i] = tr_y_yyyy_y[i] * fe_0 + tr_y_yyyyz_y[i] * pa_z[i];

        tr_y_yyyyzz_z[i] = 3.0 * tr_y_yyzz_z[i] * fe_0 + ts_yyyzz_z[i] * fe_0 + tr_y_yyyzz_z[i] * pa_y[i];
    }

    // Set up 156-159 components of targeted buffer : IP

    auto tr_y_yyyzzz_x = pbuffer.data(idx_dip_ip + 156);

    auto tr_y_yyyzzz_y = pbuffer.data(idx_dip_ip + 157);

    auto tr_y_yyyzzz_z = pbuffer.data(idx_dip_ip + 158);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_yyyz_x,   \
                             tr_y_yyyz_y,   \
                             tr_y_yyyzz_x,  \
                             tr_y_yyyzz_y,  \
                             tr_y_yyyzzz_x, \
                             tr_y_yyyzzz_y, \
                             tr_y_yyyzzz_z, \
                             tr_y_yyzzz_z,  \
                             tr_y_yzzz_z,   \
                             ts_yyzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzzz_x[i] = 2.0 * tr_y_yyyz_x[i] * fe_0 + tr_y_yyyzz_x[i] * pa_z[i];

        tr_y_yyyzzz_y[i] = 2.0 * tr_y_yyyz_y[i] * fe_0 + tr_y_yyyzz_y[i] * pa_z[i];

        tr_y_yyyzzz_z[i] = 2.0 * tr_y_yzzz_z[i] * fe_0 + ts_yyzzz_z[i] * fe_0 + tr_y_yyzzz_z[i] * pa_y[i];
    }

    // Set up 159-162 components of targeted buffer : IP

    auto tr_y_yyzzzz_x = pbuffer.data(idx_dip_ip + 159);

    auto tr_y_yyzzzz_y = pbuffer.data(idx_dip_ip + 160);

    auto tr_y_yyzzzz_z = pbuffer.data(idx_dip_ip + 161);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_yyzz_x,   \
                             tr_y_yyzz_y,   \
                             tr_y_yyzzz_x,  \
                             tr_y_yyzzz_y,  \
                             tr_y_yyzzzz_x, \
                             tr_y_yyzzzz_y, \
                             tr_y_yyzzzz_z, \
                             tr_y_yzzzz_z,  \
                             tr_y_zzzz_z,   \
                             ts_yzzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzzz_x[i] = 3.0 * tr_y_yyzz_x[i] * fe_0 + tr_y_yyzzz_x[i] * pa_z[i];

        tr_y_yyzzzz_y[i] = 3.0 * tr_y_yyzz_y[i] * fe_0 + tr_y_yyzzz_y[i] * pa_z[i];

        tr_y_yyzzzz_z[i] = tr_y_zzzz_z[i] * fe_0 + ts_yzzzz_z[i] * fe_0 + tr_y_yzzzz_z[i] * pa_y[i];
    }

    // Set up 162-165 components of targeted buffer : IP

    auto tr_y_yzzzzz_x = pbuffer.data(idx_dip_ip + 162);

    auto tr_y_yzzzzz_y = pbuffer.data(idx_dip_ip + 163);

    auto tr_y_yzzzzz_z = pbuffer.data(idx_dip_ip + 164);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_yzzz_y,   \
                             tr_y_yzzzz_y,  \
                             tr_y_yzzzzz_x, \
                             tr_y_yzzzzz_y, \
                             tr_y_yzzzzz_z, \
                             tr_y_zzzzz_x,  \
                             tr_y_zzzzz_z,  \
                             ts_zzzzz_x,    \
                             ts_zzzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzzz_x[i] = ts_zzzzz_x[i] * fe_0 + tr_y_zzzzz_x[i] * pa_y[i];

        tr_y_yzzzzz_y[i] = 4.0 * tr_y_yzzz_y[i] * fe_0 + tr_y_yzzzz_y[i] * pa_z[i];

        tr_y_yzzzzz_z[i] = ts_zzzzz_z[i] * fe_0 + tr_y_zzzzz_z[i] * pa_y[i];
    }

    // Set up 165-168 components of targeted buffer : IP

    auto tr_y_zzzzzz_x = pbuffer.data(idx_dip_ip + 165);

    auto tr_y_zzzzzz_y = pbuffer.data(idx_dip_ip + 166);

    auto tr_y_zzzzzz_z = pbuffer.data(idx_dip_ip + 167);

#pragma omp simd aligned(pa_z,              \
                             tr_y_zzzz_x,   \
                             tr_y_zzzz_y,   \
                             tr_y_zzzz_z,   \
                             tr_y_zzzzz_0,  \
                             tr_y_zzzzz_x,  \
                             tr_y_zzzzz_y,  \
                             tr_y_zzzzz_z,  \
                             tr_y_zzzzzz_x, \
                             tr_y_zzzzzz_y, \
                             tr_y_zzzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzzz_x[i] = 5.0 * tr_y_zzzz_x[i] * fe_0 + tr_y_zzzzz_x[i] * pa_z[i];

        tr_y_zzzzzz_y[i] = 5.0 * tr_y_zzzz_y[i] * fe_0 + tr_y_zzzzz_y[i] * pa_z[i];

        tr_y_zzzzzz_z[i] = 5.0 * tr_y_zzzz_z[i] * fe_0 + tr_y_zzzzz_0[i] * fe_0 + tr_y_zzzzz_z[i] * pa_z[i];
    }

    // Set up 168-171 components of targeted buffer : IP

    auto tr_z_xxxxxx_x = pbuffer.data(idx_dip_ip + 168);

    auto tr_z_xxxxxx_y = pbuffer.data(idx_dip_ip + 169);

    auto tr_z_xxxxxx_z = pbuffer.data(idx_dip_ip + 170);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xxxx_x,   \
                             tr_z_xxxx_y,   \
                             tr_z_xxxx_z,   \
                             tr_z_xxxxx_0,  \
                             tr_z_xxxxx_x,  \
                             tr_z_xxxxx_y,  \
                             tr_z_xxxxx_z,  \
                             tr_z_xxxxxx_x, \
                             tr_z_xxxxxx_y, \
                             tr_z_xxxxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxx_x[i] = 5.0 * tr_z_xxxx_x[i] * fe_0 + tr_z_xxxxx_0[i] * fe_0 + tr_z_xxxxx_x[i] * pa_x[i];

        tr_z_xxxxxx_y[i] = 5.0 * tr_z_xxxx_y[i] * fe_0 + tr_z_xxxxx_y[i] * pa_x[i];

        tr_z_xxxxxx_z[i] = 5.0 * tr_z_xxxx_z[i] * fe_0 + tr_z_xxxxx_z[i] * pa_x[i];
    }

    // Set up 171-174 components of targeted buffer : IP

    auto tr_z_xxxxxy_x = pbuffer.data(idx_dip_ip + 171);

    auto tr_z_xxxxxy_y = pbuffer.data(idx_dip_ip + 172);

    auto tr_z_xxxxxy_z = pbuffer.data(idx_dip_ip + 173);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxx_x, tr_z_xxxxx_z, tr_z_xxxxxy_x, tr_z_xxxxxy_y, tr_z_xxxxxy_z, tr_z_xxxxy_y, tr_z_xxxy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxy_x[i] = tr_z_xxxxx_x[i] * pa_y[i];

        tr_z_xxxxxy_y[i] = 4.0 * tr_z_xxxy_y[i] * fe_0 + tr_z_xxxxy_y[i] * pa_x[i];

        tr_z_xxxxxy_z[i] = tr_z_xxxxx_z[i] * pa_y[i];
    }

    // Set up 174-177 components of targeted buffer : IP

    auto tr_z_xxxxxz_x = pbuffer.data(idx_dip_ip + 174);

    auto tr_z_xxxxxz_y = pbuffer.data(idx_dip_ip + 175);

    auto tr_z_xxxxxz_z = pbuffer.data(idx_dip_ip + 176);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_z_xxxxx_x,  \
                             tr_z_xxxxxz_x, \
                             tr_z_xxxxxz_y, \
                             tr_z_xxxxxz_z, \
                             tr_z_xxxxz_y,  \
                             tr_z_xxxxz_z,  \
                             tr_z_xxxz_y,   \
                             tr_z_xxxz_z,   \
                             ts_xxxxx_x,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxxz_x[i] = ts_xxxxx_x[i] * fe_0 + tr_z_xxxxx_x[i] * pa_z[i];

        tr_z_xxxxxz_y[i] = 4.0 * tr_z_xxxz_y[i] * fe_0 + tr_z_xxxxz_y[i] * pa_x[i];

        tr_z_xxxxxz_z[i] = 4.0 * tr_z_xxxz_z[i] * fe_0 + tr_z_xxxxz_z[i] * pa_x[i];
    }

    // Set up 177-180 components of targeted buffer : IP

    auto tr_z_xxxxyy_x = pbuffer.data(idx_dip_ip + 177);

    auto tr_z_xxxxyy_y = pbuffer.data(idx_dip_ip + 178);

    auto tr_z_xxxxyy_z = pbuffer.data(idx_dip_ip + 179);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxxx_x,   \
                             tr_z_xxxxy_x,  \
                             tr_z_xxxxyy_x, \
                             tr_z_xxxxyy_y, \
                             tr_z_xxxxyy_z, \
                             tr_z_xxxyy_y,  \
                             tr_z_xxxyy_z,  \
                             tr_z_xxyy_y,   \
                             tr_z_xxyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyy_x[i] = tr_z_xxxx_x[i] * fe_0 + tr_z_xxxxy_x[i] * pa_y[i];

        tr_z_xxxxyy_y[i] = 3.0 * tr_z_xxyy_y[i] * fe_0 + tr_z_xxxyy_y[i] * pa_x[i];

        tr_z_xxxxyy_z[i] = 3.0 * tr_z_xxyy_z[i] * fe_0 + tr_z_xxxyy_z[i] * pa_x[i];
    }

    // Set up 180-183 components of targeted buffer : IP

    auto tr_z_xxxxyz_x = pbuffer.data(idx_dip_ip + 180);

    auto tr_z_xxxxyz_y = pbuffer.data(idx_dip_ip + 181);

    auto tr_z_xxxxyz_z = pbuffer.data(idx_dip_ip + 182);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxxxyz_x, tr_z_xxxxyz_y, tr_z_xxxxyz_z, tr_z_xxxxz_x, tr_z_xxxxz_z, tr_z_xxxyz_y, tr_z_xxyz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxyz_x[i] = tr_z_xxxxz_x[i] * pa_y[i];

        tr_z_xxxxyz_y[i] = 3.0 * tr_z_xxyz_y[i] * fe_0 + tr_z_xxxyz_y[i] * pa_x[i];

        tr_z_xxxxyz_z[i] = tr_z_xxxxz_z[i] * pa_y[i];
    }

    // Set up 183-186 components of targeted buffer : IP

    auto tr_z_xxxxzz_x = pbuffer.data(idx_dip_ip + 183);

    auto tr_z_xxxxzz_y = pbuffer.data(idx_dip_ip + 184);

    auto tr_z_xxxxzz_z = pbuffer.data(idx_dip_ip + 185);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xxxxzz_x, \
                             tr_z_xxxxzz_y, \
                             tr_z_xxxxzz_z, \
                             tr_z_xxxzz_0,  \
                             tr_z_xxxzz_x,  \
                             tr_z_xxxzz_y,  \
                             tr_z_xxxzz_z,  \
                             tr_z_xxzz_x,   \
                             tr_z_xxzz_y,   \
                             tr_z_xxzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxzz_x[i] = 3.0 * tr_z_xxzz_x[i] * fe_0 + tr_z_xxxzz_0[i] * fe_0 + tr_z_xxxzz_x[i] * pa_x[i];

        tr_z_xxxxzz_y[i] = 3.0 * tr_z_xxzz_y[i] * fe_0 + tr_z_xxxzz_y[i] * pa_x[i];

        tr_z_xxxxzz_z[i] = 3.0 * tr_z_xxzz_z[i] * fe_0 + tr_z_xxxzz_z[i] * pa_x[i];
    }

    // Set up 186-189 components of targeted buffer : IP

    auto tr_z_xxxyyy_x = pbuffer.data(idx_dip_ip + 186);

    auto tr_z_xxxyyy_y = pbuffer.data(idx_dip_ip + 187);

    auto tr_z_xxxyyy_z = pbuffer.data(idx_dip_ip + 188);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxxy_x,   \
                             tr_z_xxxyy_x,  \
                             tr_z_xxxyyy_x, \
                             tr_z_xxxyyy_y, \
                             tr_z_xxxyyy_z, \
                             tr_z_xxyyy_y,  \
                             tr_z_xxyyy_z,  \
                             tr_z_xyyy_y,   \
                             tr_z_xyyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyy_x[i] = 2.0 * tr_z_xxxy_x[i] * fe_0 + tr_z_xxxyy_x[i] * pa_y[i];

        tr_z_xxxyyy_y[i] = 2.0 * tr_z_xyyy_y[i] * fe_0 + tr_z_xxyyy_y[i] * pa_x[i];

        tr_z_xxxyyy_z[i] = 2.0 * tr_z_xyyy_z[i] * fe_0 + tr_z_xxyyy_z[i] * pa_x[i];
    }

    // Set up 189-192 components of targeted buffer : IP

    auto tr_z_xxxyyz_x = pbuffer.data(idx_dip_ip + 189);

    auto tr_z_xxxyyz_y = pbuffer.data(idx_dip_ip + 190);

    auto tr_z_xxxyyz_z = pbuffer.data(idx_dip_ip + 191);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxxyyz_x, \
                             tr_z_xxxyyz_y, \
                             tr_z_xxxyyz_z, \
                             tr_z_xxxyz_x,  \
                             tr_z_xxxz_x,   \
                             tr_z_xxyyz_y,  \
                             tr_z_xxyyz_z,  \
                             tr_z_xyyz_y,   \
                             tr_z_xyyz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyyz_x[i] = tr_z_xxxz_x[i] * fe_0 + tr_z_xxxyz_x[i] * pa_y[i];

        tr_z_xxxyyz_y[i] = 2.0 * tr_z_xyyz_y[i] * fe_0 + tr_z_xxyyz_y[i] * pa_x[i];

        tr_z_xxxyyz_z[i] = 2.0 * tr_z_xyyz_z[i] * fe_0 + tr_z_xxyyz_z[i] * pa_x[i];
    }

    // Set up 192-195 components of targeted buffer : IP

    auto tr_z_xxxyzz_x = pbuffer.data(idx_dip_ip + 192);

    auto tr_z_xxxyzz_y = pbuffer.data(idx_dip_ip + 193);

    auto tr_z_xxxyzz_z = pbuffer.data(idx_dip_ip + 194);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxxyzz_x, tr_z_xxxyzz_y, tr_z_xxxyzz_z, tr_z_xxxzz_x, tr_z_xxxzz_z, tr_z_xxyzz_y, tr_z_xyzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyzz_x[i] = tr_z_xxxzz_x[i] * pa_y[i];

        tr_z_xxxyzz_y[i] = 2.0 * tr_z_xyzz_y[i] * fe_0 + tr_z_xxyzz_y[i] * pa_x[i];

        tr_z_xxxyzz_z[i] = tr_z_xxxzz_z[i] * pa_y[i];
    }

    // Set up 195-198 components of targeted buffer : IP

    auto tr_z_xxxzzz_x = pbuffer.data(idx_dip_ip + 195);

    auto tr_z_xxxzzz_y = pbuffer.data(idx_dip_ip + 196);

    auto tr_z_xxxzzz_z = pbuffer.data(idx_dip_ip + 197);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xxxzzz_x, \
                             tr_z_xxxzzz_y, \
                             tr_z_xxxzzz_z, \
                             tr_z_xxzzz_0,  \
                             tr_z_xxzzz_x,  \
                             tr_z_xxzzz_y,  \
                             tr_z_xxzzz_z,  \
                             tr_z_xzzz_x,   \
                             tr_z_xzzz_y,   \
                             tr_z_xzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzzz_x[i] = 2.0 * tr_z_xzzz_x[i] * fe_0 + tr_z_xxzzz_0[i] * fe_0 + tr_z_xxzzz_x[i] * pa_x[i];

        tr_z_xxxzzz_y[i] = 2.0 * tr_z_xzzz_y[i] * fe_0 + tr_z_xxzzz_y[i] * pa_x[i];

        tr_z_xxxzzz_z[i] = 2.0 * tr_z_xzzz_z[i] * fe_0 + tr_z_xxzzz_z[i] * pa_x[i];
    }

    // Set up 198-201 components of targeted buffer : IP

    auto tr_z_xxyyyy_x = pbuffer.data(idx_dip_ip + 198);

    auto tr_z_xxyyyy_y = pbuffer.data(idx_dip_ip + 199);

    auto tr_z_xxyyyy_z = pbuffer.data(idx_dip_ip + 200);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxyy_x,   \
                             tr_z_xxyyy_x,  \
                             tr_z_xxyyyy_x, \
                             tr_z_xxyyyy_y, \
                             tr_z_xxyyyy_z, \
                             tr_z_xyyyy_y,  \
                             tr_z_xyyyy_z,  \
                             tr_z_yyyy_y,   \
                             tr_z_yyyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyy_x[i] = 3.0 * tr_z_xxyy_x[i] * fe_0 + tr_z_xxyyy_x[i] * pa_y[i];

        tr_z_xxyyyy_y[i] = tr_z_yyyy_y[i] * fe_0 + tr_z_xyyyy_y[i] * pa_x[i];

        tr_z_xxyyyy_z[i] = tr_z_yyyy_z[i] * fe_0 + tr_z_xyyyy_z[i] * pa_x[i];
    }

    // Set up 201-204 components of targeted buffer : IP

    auto tr_z_xxyyyz_x = pbuffer.data(idx_dip_ip + 201);

    auto tr_z_xxyyyz_y = pbuffer.data(idx_dip_ip + 202);

    auto tr_z_xxyyyz_z = pbuffer.data(idx_dip_ip + 203);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxyyyz_x, \
                             tr_z_xxyyyz_y, \
                             tr_z_xxyyyz_z, \
                             tr_z_xxyyz_x,  \
                             tr_z_xxyz_x,   \
                             tr_z_xyyyz_y,  \
                             tr_z_xyyyz_z,  \
                             tr_z_yyyz_y,   \
                             tr_z_yyyz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyyz_x[i] = 2.0 * tr_z_xxyz_x[i] * fe_0 + tr_z_xxyyz_x[i] * pa_y[i];

        tr_z_xxyyyz_y[i] = tr_z_yyyz_y[i] * fe_0 + tr_z_xyyyz_y[i] * pa_x[i];

        tr_z_xxyyyz_z[i] = tr_z_yyyz_z[i] * fe_0 + tr_z_xyyyz_z[i] * pa_x[i];
    }

    // Set up 204-207 components of targeted buffer : IP

    auto tr_z_xxyyzz_x = pbuffer.data(idx_dip_ip + 204);

    auto tr_z_xxyyzz_y = pbuffer.data(idx_dip_ip + 205);

    auto tr_z_xxyyzz_z = pbuffer.data(idx_dip_ip + 206);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxyyzz_x, \
                             tr_z_xxyyzz_y, \
                             tr_z_xxyyzz_z, \
                             tr_z_xxyzz_x,  \
                             tr_z_xxzz_x,   \
                             tr_z_xyyzz_y,  \
                             tr_z_xyyzz_z,  \
                             tr_z_yyzz_y,   \
                             tr_z_yyzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyzz_x[i] = tr_z_xxzz_x[i] * fe_0 + tr_z_xxyzz_x[i] * pa_y[i];

        tr_z_xxyyzz_y[i] = tr_z_yyzz_y[i] * fe_0 + tr_z_xyyzz_y[i] * pa_x[i];

        tr_z_xxyyzz_z[i] = tr_z_yyzz_z[i] * fe_0 + tr_z_xyyzz_z[i] * pa_x[i];
    }

    // Set up 207-210 components of targeted buffer : IP

    auto tr_z_xxyzzz_x = pbuffer.data(idx_dip_ip + 207);

    auto tr_z_xxyzzz_y = pbuffer.data(idx_dip_ip + 208);

    auto tr_z_xxyzzz_z = pbuffer.data(idx_dip_ip + 209);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxyzzz_x, tr_z_xxyzzz_y, tr_z_xxyzzz_z, tr_z_xxzzz_x, tr_z_xxzzz_z, tr_z_xyzzz_y, tr_z_yzzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzzz_x[i] = tr_z_xxzzz_x[i] * pa_y[i];

        tr_z_xxyzzz_y[i] = tr_z_yzzz_y[i] * fe_0 + tr_z_xyzzz_y[i] * pa_x[i];

        tr_z_xxyzzz_z[i] = tr_z_xxzzz_z[i] * pa_y[i];
    }

    // Set up 210-213 components of targeted buffer : IP

    auto tr_z_xxzzzz_x = pbuffer.data(idx_dip_ip + 210);

    auto tr_z_xxzzzz_y = pbuffer.data(idx_dip_ip + 211);

    auto tr_z_xxzzzz_z = pbuffer.data(idx_dip_ip + 212);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xxzzzz_x, \
                             tr_z_xxzzzz_y, \
                             tr_z_xxzzzz_z, \
                             tr_z_xzzzz_0,  \
                             tr_z_xzzzz_x,  \
                             tr_z_xzzzz_y,  \
                             tr_z_xzzzz_z,  \
                             tr_z_zzzz_x,   \
                             tr_z_zzzz_y,   \
                             tr_z_zzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzzz_x[i] = tr_z_zzzz_x[i] * fe_0 + tr_z_xzzzz_0[i] * fe_0 + tr_z_xzzzz_x[i] * pa_x[i];

        tr_z_xxzzzz_y[i] = tr_z_zzzz_y[i] * fe_0 + tr_z_xzzzz_y[i] * pa_x[i];

        tr_z_xxzzzz_z[i] = tr_z_zzzz_z[i] * fe_0 + tr_z_xzzzz_z[i] * pa_x[i];
    }

    // Set up 213-216 components of targeted buffer : IP

    auto tr_z_xyyyyy_x = pbuffer.data(idx_dip_ip + 213);

    auto tr_z_xyyyyy_y = pbuffer.data(idx_dip_ip + 214);

    auto tr_z_xyyyyy_z = pbuffer.data(idx_dip_ip + 215);

#pragma omp simd aligned(pa_x, tr_z_xyyyyy_x, tr_z_xyyyyy_y, tr_z_xyyyyy_z, tr_z_yyyyy_0, tr_z_yyyyy_x, tr_z_yyyyy_y, tr_z_yyyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyy_x[i] = tr_z_yyyyy_0[i] * fe_0 + tr_z_yyyyy_x[i] * pa_x[i];

        tr_z_xyyyyy_y[i] = tr_z_yyyyy_y[i] * pa_x[i];

        tr_z_xyyyyy_z[i] = tr_z_yyyyy_z[i] * pa_x[i];
    }

    // Set up 216-219 components of targeted buffer : IP

    auto tr_z_xyyyyz_x = pbuffer.data(idx_dip_ip + 216);

    auto tr_z_xyyyyz_y = pbuffer.data(idx_dip_ip + 217);

    auto tr_z_xyyyyz_z = pbuffer.data(idx_dip_ip + 218);

#pragma omp simd aligned(pa_x, tr_z_xyyyyz_x, tr_z_xyyyyz_y, tr_z_xyyyyz_z, tr_z_yyyyz_0, tr_z_yyyyz_x, tr_z_yyyyz_y, tr_z_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyyz_x[i] = tr_z_yyyyz_0[i] * fe_0 + tr_z_yyyyz_x[i] * pa_x[i];

        tr_z_xyyyyz_y[i] = tr_z_yyyyz_y[i] * pa_x[i];

        tr_z_xyyyyz_z[i] = tr_z_yyyyz_z[i] * pa_x[i];
    }

    // Set up 219-222 components of targeted buffer : IP

    auto tr_z_xyyyzz_x = pbuffer.data(idx_dip_ip + 219);

    auto tr_z_xyyyzz_y = pbuffer.data(idx_dip_ip + 220);

    auto tr_z_xyyyzz_z = pbuffer.data(idx_dip_ip + 221);

#pragma omp simd aligned(pa_x, tr_z_xyyyzz_x, tr_z_xyyyzz_y, tr_z_xyyyzz_z, tr_z_yyyzz_0, tr_z_yyyzz_x, tr_z_yyyzz_y, tr_z_yyyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyzz_x[i] = tr_z_yyyzz_0[i] * fe_0 + tr_z_yyyzz_x[i] * pa_x[i];

        tr_z_xyyyzz_y[i] = tr_z_yyyzz_y[i] * pa_x[i];

        tr_z_xyyyzz_z[i] = tr_z_yyyzz_z[i] * pa_x[i];
    }

    // Set up 222-225 components of targeted buffer : IP

    auto tr_z_xyyzzz_x = pbuffer.data(idx_dip_ip + 222);

    auto tr_z_xyyzzz_y = pbuffer.data(idx_dip_ip + 223);

    auto tr_z_xyyzzz_z = pbuffer.data(idx_dip_ip + 224);

#pragma omp simd aligned(pa_x, tr_z_xyyzzz_x, tr_z_xyyzzz_y, tr_z_xyyzzz_z, tr_z_yyzzz_0, tr_z_yyzzz_x, tr_z_yyzzz_y, tr_z_yyzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzzz_x[i] = tr_z_yyzzz_0[i] * fe_0 + tr_z_yyzzz_x[i] * pa_x[i];

        tr_z_xyyzzz_y[i] = tr_z_yyzzz_y[i] * pa_x[i];

        tr_z_xyyzzz_z[i] = tr_z_yyzzz_z[i] * pa_x[i];
    }

    // Set up 225-228 components of targeted buffer : IP

    auto tr_z_xyzzzz_x = pbuffer.data(idx_dip_ip + 225);

    auto tr_z_xyzzzz_y = pbuffer.data(idx_dip_ip + 226);

    auto tr_z_xyzzzz_z = pbuffer.data(idx_dip_ip + 227);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xyzzzz_x, tr_z_xyzzzz_y, tr_z_xyzzzz_z, tr_z_xzzzz_x, tr_z_yzzzz_y, tr_z_yzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_z_xyzzzz_x[i] = tr_z_xzzzz_x[i] * pa_y[i];

        tr_z_xyzzzz_y[i] = tr_z_yzzzz_y[i] * pa_x[i];

        tr_z_xyzzzz_z[i] = tr_z_yzzzz_z[i] * pa_x[i];
    }

    // Set up 228-231 components of targeted buffer : IP

    auto tr_z_xzzzzz_x = pbuffer.data(idx_dip_ip + 228);

    auto tr_z_xzzzzz_y = pbuffer.data(idx_dip_ip + 229);

    auto tr_z_xzzzzz_z = pbuffer.data(idx_dip_ip + 230);

#pragma omp simd aligned(pa_x, tr_z_xzzzzz_x, tr_z_xzzzzz_y, tr_z_xzzzzz_z, tr_z_zzzzz_0, tr_z_zzzzz_x, tr_z_zzzzz_y, tr_z_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzzz_x[i] = tr_z_zzzzz_0[i] * fe_0 + tr_z_zzzzz_x[i] * pa_x[i];

        tr_z_xzzzzz_y[i] = tr_z_zzzzz_y[i] * pa_x[i];

        tr_z_xzzzzz_z[i] = tr_z_zzzzz_z[i] * pa_x[i];
    }

    // Set up 231-234 components of targeted buffer : IP

    auto tr_z_yyyyyy_x = pbuffer.data(idx_dip_ip + 231);

    auto tr_z_yyyyyy_y = pbuffer.data(idx_dip_ip + 232);

    auto tr_z_yyyyyy_z = pbuffer.data(idx_dip_ip + 233);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yyyy_x,   \
                             tr_z_yyyy_y,   \
                             tr_z_yyyy_z,   \
                             tr_z_yyyyy_0,  \
                             tr_z_yyyyy_x,  \
                             tr_z_yyyyy_y,  \
                             tr_z_yyyyy_z,  \
                             tr_z_yyyyyy_x, \
                             tr_z_yyyyyy_y, \
                             tr_z_yyyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyy_x[i] = 5.0 * tr_z_yyyy_x[i] * fe_0 + tr_z_yyyyy_x[i] * pa_y[i];

        tr_z_yyyyyy_y[i] = 5.0 * tr_z_yyyy_y[i] * fe_0 + tr_z_yyyyy_0[i] * fe_0 + tr_z_yyyyy_y[i] * pa_y[i];

        tr_z_yyyyyy_z[i] = 5.0 * tr_z_yyyy_z[i] * fe_0 + tr_z_yyyyy_z[i] * pa_y[i];
    }

    // Set up 234-237 components of targeted buffer : IP

    auto tr_z_yyyyyz_x = pbuffer.data(idx_dip_ip + 234);

    auto tr_z_yyyyyz_y = pbuffer.data(idx_dip_ip + 235);

    auto tr_z_yyyyyz_z = pbuffer.data(idx_dip_ip + 236);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_z_yyyyy_y,  \
                             tr_z_yyyyyz_x, \
                             tr_z_yyyyyz_y, \
                             tr_z_yyyyyz_z, \
                             tr_z_yyyyz_x,  \
                             tr_z_yyyyz_z,  \
                             tr_z_yyyz_x,   \
                             tr_z_yyyz_z,   \
                             ts_yyyyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyyz_x[i] = 4.0 * tr_z_yyyz_x[i] * fe_0 + tr_z_yyyyz_x[i] * pa_y[i];

        tr_z_yyyyyz_y[i] = ts_yyyyy_y[i] * fe_0 + tr_z_yyyyy_y[i] * pa_z[i];

        tr_z_yyyyyz_z[i] = 4.0 * tr_z_yyyz_z[i] * fe_0 + tr_z_yyyyz_z[i] * pa_y[i];
    }

    // Set up 237-240 components of targeted buffer : IP

    auto tr_z_yyyyzz_x = pbuffer.data(idx_dip_ip + 237);

    auto tr_z_yyyyzz_y = pbuffer.data(idx_dip_ip + 238);

    auto tr_z_yyyyzz_z = pbuffer.data(idx_dip_ip + 239);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yyyyzz_x, \
                             tr_z_yyyyzz_y, \
                             tr_z_yyyyzz_z, \
                             tr_z_yyyzz_0,  \
                             tr_z_yyyzz_x,  \
                             tr_z_yyyzz_y,  \
                             tr_z_yyyzz_z,  \
                             tr_z_yyzz_x,   \
                             tr_z_yyzz_y,   \
                             tr_z_yyzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyzz_x[i] = 3.0 * tr_z_yyzz_x[i] * fe_0 + tr_z_yyyzz_x[i] * pa_y[i];

        tr_z_yyyyzz_y[i] = 3.0 * tr_z_yyzz_y[i] * fe_0 + tr_z_yyyzz_0[i] * fe_0 + tr_z_yyyzz_y[i] * pa_y[i];

        tr_z_yyyyzz_z[i] = 3.0 * tr_z_yyzz_z[i] * fe_0 + tr_z_yyyzz_z[i] * pa_y[i];
    }

    // Set up 240-243 components of targeted buffer : IP

    auto tr_z_yyyzzz_x = pbuffer.data(idx_dip_ip + 240);

    auto tr_z_yyyzzz_y = pbuffer.data(idx_dip_ip + 241);

    auto tr_z_yyyzzz_z = pbuffer.data(idx_dip_ip + 242);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yyyzzz_x, \
                             tr_z_yyyzzz_y, \
                             tr_z_yyyzzz_z, \
                             tr_z_yyzzz_0,  \
                             tr_z_yyzzz_x,  \
                             tr_z_yyzzz_y,  \
                             tr_z_yyzzz_z,  \
                             tr_z_yzzz_x,   \
                             tr_z_yzzz_y,   \
                             tr_z_yzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzzz_x[i] = 2.0 * tr_z_yzzz_x[i] * fe_0 + tr_z_yyzzz_x[i] * pa_y[i];

        tr_z_yyyzzz_y[i] = 2.0 * tr_z_yzzz_y[i] * fe_0 + tr_z_yyzzz_0[i] * fe_0 + tr_z_yyzzz_y[i] * pa_y[i];

        tr_z_yyyzzz_z[i] = 2.0 * tr_z_yzzz_z[i] * fe_0 + tr_z_yyzzz_z[i] * pa_y[i];
    }

    // Set up 243-246 components of targeted buffer : IP

    auto tr_z_yyzzzz_x = pbuffer.data(idx_dip_ip + 243);

    auto tr_z_yyzzzz_y = pbuffer.data(idx_dip_ip + 244);

    auto tr_z_yyzzzz_z = pbuffer.data(idx_dip_ip + 245);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yyzzzz_x, \
                             tr_z_yyzzzz_y, \
                             tr_z_yyzzzz_z, \
                             tr_z_yzzzz_0,  \
                             tr_z_yzzzz_x,  \
                             tr_z_yzzzz_y,  \
                             tr_z_yzzzz_z,  \
                             tr_z_zzzz_x,   \
                             tr_z_zzzz_y,   \
                             tr_z_zzzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzzz_x[i] = tr_z_zzzz_x[i] * fe_0 + tr_z_yzzzz_x[i] * pa_y[i];

        tr_z_yyzzzz_y[i] = tr_z_zzzz_y[i] * fe_0 + tr_z_yzzzz_0[i] * fe_0 + tr_z_yzzzz_y[i] * pa_y[i];

        tr_z_yyzzzz_z[i] = tr_z_zzzz_z[i] * fe_0 + tr_z_yzzzz_z[i] * pa_y[i];
    }

    // Set up 246-249 components of targeted buffer : IP

    auto tr_z_yzzzzz_x = pbuffer.data(idx_dip_ip + 246);

    auto tr_z_yzzzzz_y = pbuffer.data(idx_dip_ip + 247);

    auto tr_z_yzzzzz_z = pbuffer.data(idx_dip_ip + 248);

#pragma omp simd aligned(pa_y, tr_z_yzzzzz_x, tr_z_yzzzzz_y, tr_z_yzzzzz_z, tr_z_zzzzz_0, tr_z_zzzzz_x, tr_z_zzzzz_y, tr_z_zzzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzzz_x[i] = tr_z_zzzzz_x[i] * pa_y[i];

        tr_z_yzzzzz_y[i] = tr_z_zzzzz_0[i] * fe_0 + tr_z_zzzzz_y[i] * pa_y[i];

        tr_z_yzzzzz_z[i] = tr_z_zzzzz_z[i] * pa_y[i];
    }

    // Set up 249-252 components of targeted buffer : IP

    auto tr_z_zzzzzz_x = pbuffer.data(idx_dip_ip + 249);

    auto tr_z_zzzzzz_y = pbuffer.data(idx_dip_ip + 250);

    auto tr_z_zzzzzz_z = pbuffer.data(idx_dip_ip + 251);

#pragma omp simd aligned(pa_z,              \
                             tr_z_zzzz_x,   \
                             tr_z_zzzz_y,   \
                             tr_z_zzzz_z,   \
                             tr_z_zzzzz_0,  \
                             tr_z_zzzzz_x,  \
                             tr_z_zzzzz_y,  \
                             tr_z_zzzzz_z,  \
                             tr_z_zzzzzz_x, \
                             tr_z_zzzzzz_y, \
                             tr_z_zzzzzz_z, \
                             ts_zzzzz_x,    \
                             ts_zzzzz_y,    \
                             ts_zzzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzzz_x[i] = 5.0 * tr_z_zzzz_x[i] * fe_0 + ts_zzzzz_x[i] * fe_0 + tr_z_zzzzz_x[i] * pa_z[i];

        tr_z_zzzzzz_y[i] = 5.0 * tr_z_zzzz_y[i] * fe_0 + ts_zzzzz_y[i] * fe_0 + tr_z_zzzzz_y[i] * pa_z[i];

        tr_z_zzzzzz_z[i] = 5.0 * tr_z_zzzz_z[i] * fe_0 + tr_z_zzzzz_0[i] * fe_0 + ts_zzzzz_z[i] * fe_0 + tr_z_zzzzz_z[i] * pa_z[i];
    }
}

}  // namespace diprec
